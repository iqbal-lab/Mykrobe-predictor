#!/usr/bin/env python
# encoding: utf-8
"""
compare.py
Created on: 2016-03-30
Author: Phelim Bradley

Compares the JSON output of mykrobe predictor

"""
import sys
import json
sys.path.append("/home/phelimb/git/Mykrobe-predictor")
import argparse
from mykrobe.utils import load_json
from mykrobe.utils import unique
from mykrobe.predict import MykrobePredictorSusceptibilityResult
parser = argparse.ArgumentParser(description='Compares the JSON output of mykrobe predictor')
parser.add_argument('truth', metavar='truth', type=str, help='truth')
parser.add_argument('--analysis', default='summary', type = str, choices = ["summary", "table", "diff"], help = "report format (options, summary, table, diff) default : summary")
parser.add_argument('--ana1', type=str, help='analyses', nargs='+', default = [])
parser.add_argument('--ana2', type=str, help='analyses', nargs='+', default = [])
parser.add_argument('--format', default='short', type = str, choices = ["short", "long"], help = "report format (options, short, long) default : short")
args = parser.parse_args()

class Stats(object):


    def __init__(self, count_comparision):
        self.count_comparision = count_comparision

    @property
    def num_samples(self):
        return self.P + self.N

    @property
    def total(self):
        return self.num_samples

    @property
    def P(self):
        return self.TP + self.FN + self.IR

    @property
    def N(self):
        return self.TN + self.FP + self.IS

    @property
    def IR(self):
        """Pheno is R predict is inconclusive """
        return self.count_comparision.get('IR', 0)

    @property
    def IS(self):
        """Pheno is S predict is inconclusive """
        return self.count_comparision.get('IS', 0)

    @property
    def FP(self):
        return self.count_comparision.get('FP', 0)

    @property
    def TP(self):
        return self.count_comparision.get('TP', 0)

    @property
    def FN(self):
        return self.count_comparision.get('FN', 0)

    @property
    def TN(self):
        return self.count_comparision.get('TN', 0)

    @property
    def unknown(self):
        return self.count_comparision.get('UNKNOWN', 0)

    @property
    def VME_str(self):
        return "%s%% (%s%%-%s%%)" % (self.VME, self.VME_LB, self.VME_UB)

    @property
    def ME_str(self):
        return "%s%% (%s%%-%s%%)" % (self.ME, self.ME_LB, self.ME_UB)

    @property
    def sensitivity_str(self):
        return "%s%% (%s%%-%s%%)" % (self.sensitivity,
                                     self.sensitivity_LB, self.sensitivity_UB)

    @property
    def specificity_str(self):
        return "%s%% (%s%%-%s%%)" % (self.specificity,
                                     self.specificity_LB, self.specificity_UB)

    @property
    def FN_str(self):
        return "%s (%s)" % (self.FN, self.P)

    @property
    def FP_str(self):
        return "%s (%s)" % (self.FP, self.N)

    def percentage(self, num, denom):
        try:
            return round(100 * (float(num) / denom), 1)
        except ZeroDivisionError:
            return -1

    @property
    def VME(self):
        return self.percentage(self.FN, self.P)

    @property
    def ME(self):
        return self.percentage(self.FP, self.N)

    @property
    def sensitivity(self):
        return self.percentage(self.TP, self.P)

    @property
    def specificity(self):
        return self.percentage(self.TN, self.N)

    @property
    def VME_UB(self):
        return round(100 * (self.VME_conf[1]), 1)

    @property
    def VME_LB(self):
        return round(100 * (self.VME_conf[0]), 1)

    @property
    def ME_UB(self):
        return round(100 * (self.ME_conf[1]), 1)

    @property
    def ME_LB(self):
        return round(100 * (self.ME_conf[0]), 1)

    @property
    def sensitivity_UB(self):
        return round(100 * (self.sensitivity_conf[1]), 1)

    @property
    def sensitivity_LB(self):
        return round(100 * (self.sensitivity_conf[0]), 1)

    @property
    def specificity_UB(self):
        return round(100 * (self.specificity_conf[1]), 1)

    @property
    def specificity_LB(self):
        return round(100 * (self.specificity_conf[0]), 1)

    @property
    def PPV(self):
        return self.percentage(self.TP , (self.TP + self.FP))

    @property
    def PPV_conf(self):
        return self.binom_interval(success=self.TP, total=(self.TP + self.FP) )      

    @property
    def PPV_UB(self):
        return round(100 * (self.PPV_conf[1]), 1)

    @property
    def PPV_LB(self):
        return round(100 * (self.PPV_conf[0]), 1) 

    @property
    def PPV_str(self):
        return "%s%% (%s%%-%s%%)" % (self.PPV, self.PPV_LB, self.PPV_UB)                         

    @property
    def NPV(self):
        return self.percentage(self.TN , (self.TN + self.FN)  )

    @property
    def NPV_conf(self):
        return self.binom_interval(success=self.TN, total=(self.TN + self.FN) )     

    @property
    def NPV_UB(self):
        return round(100 * (self.NPV_conf[1]), 1)

    @property
    def NPV_LB(self):
        return round(100 * (self.NPV_conf[0]), 1)   

    @property
    def NPV_str(self):
        return "%s%% (%s%%-%s%%)" % (self.NPV, self.NPV_LB, self.NPV_UB)

    @property
    def VME_conf(self):
        return self.binom_interval(success=self.FN, total=self.P)

    @property
    def ME_conf(self):
        return self.binom_interval(success=self.FP, total=self.N)

    @property
    def sensitivity_conf(self):
        return self.binom_interval(success=self.TP, total=self.P)

    @property
    def specificity_conf(self):
        return self.binom_interval(success=self.TN, total=self.N)

    def binom_interval(self, success, total, confint=0.95):
        quantile = (1 - confint) / 2.
        # lower = scipy.stats.beta.ppf(quantile, success, total - success + 1)
        # upper = scipy.stats.beta.ppf(1 - quantile, success + 1, total - success)
        # return (lower, upper)
        return (0, 1)

    @property
    def row_long_header(self):
        header = ["Total", "TP", "FP", "P", "TN", "FN", "N", "VME", "ME", "Sens", "Spec",
                  "VME_LB", "VME_UB", "ME_LB", "ME_UB", "Sens_LB", "Sens_UB", "Spec_LB", "Spec_UB",
                    "PPV", "PPV_LB", "PPV_UB",
                    "NPV", "NPV_LB", "NPV_UB"
                  ]
        return header

    @property
    def row_long(self):
		return [ self.num_samples, self.TP, self.FP, self.P,
                    self.TN, self.FN, self.N, self.VME, self.ME, self.sensitivity,
                    self.specificity, self.VME_LB, self.VME_UB, self.ME_LB, self.ME_UB,
                    self.sensitivity_LB, self.sensitivity_UB, self.specificity_LB,
                    self.specificity_UB, 
                    self.PPV, self.PPV_LB, self.PPV_UB,
                    self.NPV, self.NPV_LB, self.NPV_UB, 
                    ]
    @property
    def row_short(self):
        return [self.num_samples, self.FN_str, self.FP_str,
                    self.VME_str, self.ME_str, self.sensitivity_str, self.specificity_str,
                    self.PPV_str, self.NPV_str]

    @property
    def row_short_header(self):
        header = [
            "Total",
            "FN(R)",
            "FP(S)",
            "VME",
            "ME",
            "sensitivity",
            "specificity",
            "PPV", "NPV"]
        return header    
    


## First generate a table with a column for each commit for each sample for each drug
def file_paths_to_combined_dict(l):
	ana = {}
	for f in l:
		try:
			data = load_json(f)
		except ValueError, e:
			sys.stderr.write(str(e) + " %s \n" % f)
		else:
			assert data.keys()[0] not in ana
			ana.update(data)
	return ana

def combined_dict_to_result_objects(data):
	out_dict = {}
	for k,v in data.items():
		out_dict[k] = MykrobePredictorSusceptibilityResult.from_json(json.dumps({"susceptibility" : v.get("susceptibility", {})}))
	return out_dict

def get_sample_ids(ana1, ana2):
	return unique(ana1.keys() + ana2.keys())

def create_comparision_table(sample_ids, truth, ana1, ana2):
	df = []
	for sample_id in sample_ids:
		indv_truth = truth.get(sample_id, MykrobePredictorSusceptibilityResult.create({}))
		indv_ana1 = ana1.get(sample_id, MykrobePredictorSusceptibilityResult.create({}))
		indv_ana2 = ana2.get(sample_id, MykrobePredictorSusceptibilityResult.create({}))
		drugs = unique(indv_truth.drugs + indv_ana1.drugs + indv_ana2.drugs)
		if drugs:
			for drug in drugs:
				if args.ana2:
					row= [sample_id, drug, indv_truth.susceptibility.get(drug, {"predict" : "NA"}).get("predict"), indv_ana1.susceptibility.get(drug, {"predict" : "NA"}).get("predict"), indv_ana2.susceptibility.get(drug, {"predict" : "NA"}).get("predict") ]
				else:
					row= [sample_id, drug, indv_truth.susceptibility.get(drug, {"predict" : "NA"}).get("predict"), indv_ana1.susceptibility.get(drug, {"predict" : "NA"}).get("predict") ]
				df.append(row)
	return df

def inc_count(d, k1, k2):
	if not k1 in d:
		d[k1] = {}
	if not k2 in d[k1]:
		d[k1][k2] = 0
	d[k1][k2] += 1
	return d


def update_comparision(comparison, drug, compare):
	comparison = inc_count(comparison, "all", compare)
	comparison = inc_count(comparison, drug, compare)	
	return comparison

def compare_analysis_to_truth(sample_ids, truth, ana):
	comparison = {}
	for sample_id in sample_ids:
		indv_truth = truth.get(sample_id, MykrobePredictorSusceptibilityResult.create({}))
		indv_ana = ana.get(sample_id, MykrobePredictorSusceptibilityResult.create({}))
		drugs = unique(indv_truth.drugs + indv_ana.drugs)
		if drugs:
			for drug in drugs:
				truth_drug_predict = indv_truth.susceptibility.get(drug, {"predict" : "NA"}).get("predict")
				ana_drug_predict = indv_ana.susceptibility.get(drug, {"predict" : "NA"}).get("predict")
				assert truth_drug_predict in ["R", "NA", "S"]
				if truth_drug_predict == "NA" or ana_drug_predict == "NA":
					compare = "UNKNOWN"
				elif truth_drug_predict == ana_drug_predict:
					if truth_drug_predict == "R":
						compare = "TP"
					elif truth_drug_predict == "S":
						compare = "TN"
				else:
					assert ana_drug_predict in ["R","S"]
					if truth_drug_predict == "R":
						compare = "FN"
					elif truth_drug_predict == "S":
						compare = "FP"					
				comparison = update_comparision(comparison, drug, compare)

	return comparison

def diff_stats(stats1, stats2):
	TP =  stats1.TP - stats2.TP
	FP =  stats1.FP - stats2.FP
	TN =  stats1.TN - stats2.TN
	FN =  stats1.FN - stats2.FN

	sensitivity = stats1.sensitivity - stats2.sensitivity
	specificity = stats1.specificity - stats2.specificity
	total = stats1.total - stats2.total

	return "Total: %+i, TP : %+i, FP : %+i, TN : %+i, FN : %+i, sensitivity : %+f%%, specificity : %+f%%" % (total, TP, FP, TN, FN, sensitivity, specificity)

## Load data
truth = load_json(args.truth)
ana1 = file_paths_to_combined_dict(args.ana1)
ana2 = file_paths_to_combined_dict(args.ana2)

truth_susceptibility = combined_dict_to_result_objects(truth)
ana1_susceptibility = combined_dict_to_result_objects(ana1)
ana2_susceptibility = combined_dict_to_result_objects(ana2)

sample_ids = get_sample_ids( ana1, ana2)

## Run analyses
if args.analysis == "table":
	## sample  drug  truth  ana1   ana2 
	## 1234    RIF   R     R     S 	
	df = create_comparision_table(sample_ids, truth_susceptibility, ana1_susceptibility, ana2_susceptibility)
	for row in df:
		print ("\t".join(row))
elif args.analysis == "summary":
	## Report a summary of each analysis vs. truth

	## ANA 1 
	## Drug TP FP 
	##

	## Ana 2 

	##
	##
	print "Ana1"
	count_comparision_ana1 = compare_analysis_to_truth(sample_ids, truth_susceptibility, ana1_susceptibility)
	if args.format == "short":
		print ("\t".join( ["Drug"] + [str(i) for i in Stats({}).row_short_header]))
		for k,v in count_comparision_ana1.items():
			stats = Stats(count_comparision = v)
			print ("\t".join( [k] + [str(i) for i in stats.row_short]))
		if args.ana2:
			print "Ana2"
			count_comparision_ana2 = compare_analysis_to_truth(sample_ids, truth_susceptibility, ana2_susceptibility)
			print ("\t".join( ["Drug"] + [str(i) for i in Stats({}).row_short_header]))
			for k,v in count_comparision_ana2.items():
				stats = Stats(count_comparision = v)
				print ("\t".join( [k] + [str(i) for i in stats.row_short]))


			## Diff ana 1 ana 2 summary
			print "diff ana2 - ana1"
			for k in count_comparision_ana1.keys():
				stats_2 = Stats(count_comparision = count_comparision_ana2[k])
				stats_1 = Stats(count_comparision = count_comparision_ana1[k])
				diff = diff_stats(stats_2,stats_1)
				print "\t".join([k, diff])


			##
			##				
	elif args.format == "long":
		print ("\t".join( ["Drug"] + [str(i) for i in Stats({}).row_long_header]))
		for k,v in count_comparision.items():
			stats = Stats(count_comparision = v)
			print ("\t".join( [k] + [str(i) for i in stats.row_long]))
		if args.ana2:
			print "Ana2"
			count_comparision = compare_analysis_to_truth(sample_ids, truth_susceptibility, ana2_susceptibility)
			print ("\t".join( ["Drug"] + [str(i) for i in Stats({}).row_long_header]))
			for k,v in count_comparision.items():
				stats = Stats(count_comparision = v)
				print ("\t".join( [k] + [str(i) for i in stats.row_long]))

			## Diff ana 1 ana 2 summary
			print "diff ana2 - ana1"
			for k in count_comparision_ana1.keys():
				stats_2 = Stats(count_comparision = count_comparision_ana2[k])
				stats_1 = Stats(count_comparision = count_comparision_ana1[k])
				diff = diff_stats(stats_2,stats_1)
				print "\t".join([k, diff])				

## Report the diference between ana1 and and2
elif args.analysis == "diff":
	if args.ana2:
		print "\t".join(["sample", "drug", "truth", "ana1", "ana2"])
	else:
		print "\t".join(["sample", "drug", "truth", "ana1"])
	for sample_id in sample_ids:
		indv_truth = truth_susceptibility.get(sample_id, MykrobePredictorSusceptibilityResult.create({}))
		indv_ana1 = ana1_susceptibility.get(sample_id, MykrobePredictorSusceptibilityResult.create({}))
		indv_ana2 = ana2_susceptibility.get(sample_id, MykrobePredictorSusceptibilityResult.create({}))
		if args.ana2:
			diff = indv_ana1.diff(indv_ana2)
			if diff:
				for drug, predict_diff in diff.items():
					truth_drug_predict = indv_truth.susceptibility.get(drug, {"predict" : "NA"}).get("predict")
					ana1_predict, ana2_predict = predict_diff["predict"]
					print "\t".join([sample_id, drug, truth_drug_predict, ana1_predict, ana2_predict])
		else:
			diff = indv_truth.diff(indv_ana1)
			if diff:
				for drug, predict_diff in diff.items():
					truth_drug_predict, ana1_predict = predict_diff["predict"]
					if truth_drug_predict != "NA" and ana1_predict != "NA":
						print "\t".join([sample_id, drug, truth_drug_predict, ana1_predict])			


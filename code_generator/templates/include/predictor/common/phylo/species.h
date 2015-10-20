/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 *  species.h
*/
#include "dB_graph.h"
#include "base_species.h"


{% for phylogroup in selfer.phylo_groups %}

void print_json_{{phylogroup.name}}_start();
void print_json_{{phylogroup.name}}_end();
char* get_ith_{{phylogroup.name}}_name(CovgInfo* covg_info, int i);

	typedef enum 
	 {
	 	{% for taxon in phylogroup.taxons %}
	 	{{taxon.enum}} = {{loop.index0}},
	 	{% endfor %}
    unknown{{phylogroup.name}}={{phylogroup.taxons | length}}
	   	} {{phylogroup.enum}} ;
	#define NUM_{{phylogroup.enum}} {{phylogroup.taxons | length}}
   	
  void map_{{phylogroup.name}}_enum_to_str({{phylogroup.enum}} sp, StrBuf* sbuf);
  void load_all_{{phylogroup.name}}_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir );
  void {{phylogroup.name}}_threshold(int* thresholds);
{% endfor %}

typedef struct
{
{% for phylogroup in selfer.phylo_groups %}
  CovgInfo* {{phylogroup.name}}_covg_info;
{% endfor %}
} SpeciesInfo;

void print_json_phylogenetics(SpeciesInfo* species_info);


SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last);

{% for phylogroup in selfer.phylo_groups %}
  void print_json_{{phylogroup.name}}(SpeciesInfo* species_info);
{% endfor %}


{% block extra %}
{% endblock %}
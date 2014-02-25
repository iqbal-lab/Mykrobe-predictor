/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  test_cmd_line.c
*/

#include <stdlib.h>
#include <string.h>

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "cmd_line.h"
#include "element.h"
#include "db_complex_genotyping.h"

void test_get_numbers_from_comma_sep_list()
{

  //simple list
  char list[]="1,2,3";

  int return_list[10];
  int max=10;
  int ans =  get_numbers_from_comma_sep_list(list, return_list, max);

  CU_ASSERT(return_list[0]==1);
  CU_ASSERT(return_list[1]==2);
  CU_ASSERT(return_list[2]==3);
  CU_ASSERT(ans==3);

  //list of just one number
  char list2[]="100";

  int return_list2[10];
  int max2=10;
  int ans2=get_numbers_from_comma_sep_list(list2, return_list2, max2);
  
  CU_ASSERT(ans2==1);
  printf("ans is %d \n", ans);
  CU_ASSERT(return_list2[0]==100);


  //list containing negative number
  char list3[]="-17";

  int return_list3[10];
  int max3=10;
  int ans3=get_numbers_from_comma_sep_list(list3,  return_list3, max3);
  
  CU_ASSERT(ans3==-1);


  //list containing negative number and other numbers
  char list4[]="10,1,2,3,4,-17";

  int return_list4[10];
  int max4=10;
  int ans4=get_numbers_from_comma_sep_list(list4,  return_list4, max4);
  
  CU_ASSERT(ans4==-1);

  
  //list containing too many numbers
  char list5[]="2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3";

  int return_list5[10];
  int max5=10;
  int ans5=get_numbers_from_comma_sep_list(list5,  return_list5, max5);
  
  CU_ASSERT(ans5==-1);



  //list containing a non-number
  char list6[]="10,1,2,3,6,a,2,132";

  int return_list6[10];
  int max6=10;
  int ans6=get_numbers_from_comma_sep_list(list6,  return_list6, max6);
  
  CU_ASSERT(ans6==-1);


  //list which is just -1 - this means ALL colours
  char list7[]="-1";
  int return_list7[NUMBER_OF_COLOURS];
  int max7 = NUMBER_OF_COLOURS;
  int ans7=get_numbers_from_comma_sep_list(list7,  return_list7, max7);
  int i;

  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      CU_ASSERT(return_list7[i]==i);
    }


  //list which is just one number with a * on the front, means ignore this colour, but take all others
  char list8[]="*0";
  int return_list8[NUMBER_OF_COLOURS];
  int max8 = NUMBER_OF_COLOURS;
  int ans8=get_numbers_from_comma_sep_list(list8,  return_list8, max8);

  CU_ASSERT(ans8==NUMBER_OF_COLOURS-1);
  for (i=1; i<NUMBER_OF_COLOURS; i++)
    {
      CU_ASSERT(return_list8[i-1]==i);
    }




  if (NUMBER_OF_COLOURS<10)
    {
      die("We need NUMBER_OF_COLOURS>=10 to run all these cmdline tests - recompile please!\n");
    }
  //list which is *1,3,4
  char list9[]="*1,3,4";
  int return_list9[NUMBER_OF_COLOURS];
  int max9 = NUMBER_OF_COLOURS;
  int ans9=get_numbers_from_comma_sep_list(list9,  return_list9, max9);

  CU_ASSERT(ans9==NUMBER_OF_COLOURS-3);
  CU_ASSERT(return_list9[0]==0);
  CU_ASSERT(return_list9[1]==2);

  for (i=2; i<NUMBER_OF_COLOURS-3; i++)
    {
      CU_ASSERT(return_list9[i]==5+i-2);
    }
  

}



void test_parse_colourinfo_argument()
{

  CmdLine cmd_line;

  char arg[]="1,2/3,4";
  int len_arg = strlen(arg);
  char text[]="dummy text";
  int which=1;
  parse_colourinfo_argument(&cmd_line, arg, len_arg, text, which);

  //printf("Num us  %d, ",  cmd_line.num_colours_in_detect_bubbles1_first_colour_list);
  CU_ASSERT(cmd_line.num_colours_in_detect_bubbles1_first_colour_list ==2);
  CU_ASSERT(cmd_line.detect_bubbles1_first_colour_list[0]==1);
  CU_ASSERT(cmd_line.detect_bubbles1_first_colour_list[1]==2);

  CU_ASSERT(cmd_line.num_colours_in_detect_bubbles1_second_colour_list==2);
  CU_ASSERT(cmd_line.detect_bubbles1_second_colour_list[0]==3);
  CU_ASSERT(cmd_line.detect_bubbles1_second_colour_list[1]==4);


  CmdLine cmd_line2;
  char arg2[]="1/234,17";
  int len_arg2=strlen(arg2);
  which=2;
  
  parse_colourinfo_argument(&cmd_line2, arg2, len_arg2, text, which);
  

  CU_ASSERT(cmd_line2.num_colours_in_detect_bubbles2_first_colour_list ==1);
  CU_ASSERT(cmd_line2.detect_bubbles2_first_colour_list[0]==1);

  CU_ASSERT(cmd_line2.num_colours_in_detect_bubbles2_second_colour_list==2);
  CU_ASSERT(cmd_line2.detect_bubbles2_second_colour_list[0]==234);
  CU_ASSERT(cmd_line2.detect_bubbles2_second_colour_list[1]==17);



  CmdLine cmd_line3;
  char arg3[]="1,28,19,278,8888/234,17,8";
  int len_arg3=strlen(arg3);
  which=2;
  parse_colourinfo_argument(&cmd_line3, arg3, len_arg3, text, which);


  CU_ASSERT(cmd_line3.num_colours_in_detect_bubbles2_first_colour_list ==5);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[0]==1);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[1]==28);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[2]==19);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[3]==278);
  CU_ASSERT(cmd_line3.detect_bubbles2_first_colour_list[4]==8888);

  CU_ASSERT(cmd_line3.num_colours_in_detect_bubbles2_second_colour_list==3);
  CU_ASSERT(cmd_line3.detect_bubbles2_second_colour_list[0]==234);
  CU_ASSERT(cmd_line3.detect_bubbles2_second_colour_list[1]==17);
  CU_ASSERT(cmd_line3.detect_bubbles2_second_colour_list[2]==8);







  CmdLine cmd_line4;
  char arg4[]="-1/*2,3";
  int len_arg4=strlen(arg4);
  which=2;
  parse_colourinfo_argument(&cmd_line4, arg4, len_arg4, text, which);


  CU_ASSERT(cmd_line4.num_colours_in_detect_bubbles2_first_colour_list ==NUMBER_OF_COLOURS);
  int i;
  for (i=0; i<NUMBER_OF_COLOURS; i++)
    {
      CU_ASSERT(cmd_line4.detect_bubbles2_first_colour_list[i]==i);
    }

  CU_ASSERT(cmd_line4.num_colours_in_detect_bubbles2_second_colour_list==NUMBER_OF_COLOURS-2);
  
  CU_ASSERT(cmd_line4.detect_bubbles2_second_colour_list[0]==0);
  CU_ASSERT(cmd_line4.detect_bubbles2_second_colour_list[1]==1);
  for (i=2; i<NUMBER_OF_COLOURS-2; i++)
    {
      CU_ASSERT(cmd_line4.detect_bubbles2_second_colour_list[i]==4+i-2);
    }


}


void test_parse_commasep_list()
{

  CmdLine cmd_line;

  char arg[]="1,2";
  int len_arg = strlen(arg);
  char text[]="dummy text";
  int which=1;
  int ret = parse_commasep_list(&cmd_line, arg, len_arg, text);

  CU_ASSERT(cmd_line.num_colours_in_pd_colour_list==2);
  CU_ASSERT(cmd_line.pd_colour_list[0]==1);
  CU_ASSERT(cmd_line.pd_colour_list[1]==2);
  CU_ASSERT(ret==0);

  CmdLine cmd_line2;
  char arg2[]="234";
  int len_arg2=strlen(arg2);
  
  ret = parse_commasep_list(&cmd_line2, arg2, len_arg2, text);

  CU_ASSERT(cmd_line2.num_colours_in_pd_colour_list==1);
  CU_ASSERT(cmd_line2.pd_colour_list[0]==234);
  CU_ASSERT(ret==0);
  


  CmdLine cmd_line3;
  char arg3[]="1,28,19,278,8888";
  int len_arg3=strlen(arg3);
  ret=parse_commasep_list(&cmd_line3, arg3, len_arg3, text);

  CU_ASSERT(cmd_line3.num_colours_in_pd_colour_list==5);
  CU_ASSERT(cmd_line3.pd_colour_list[0]==1);
  CU_ASSERT(cmd_line3.pd_colour_list[1]==28);
  CU_ASSERT(cmd_line3.pd_colour_list[2]==19);
  CU_ASSERT(cmd_line3.pd_colour_list[3]==278);
  CU_ASSERT(cmd_line3.pd_colour_list[4]==8888);
  CU_ASSERT(ret==0);


}



void test_parse_genotype_site_argument()
{
  
  char arg[]="1,2[3[3[2,4[zam/iqbal[UNCLEANED"; //means gt colours 1 and 2, and that 3 is the ref-minus-site colours, and there are 3 alleles, and that we must gt genotypes 2 to 4 of the 3 choose 2  = 6 possible genotypes.
  int len_arg = strlen(arg);
  int list_cols[NUMBER_OF_COLOURS];
  int num_cols_to_gt=0;
  int ref_min_site_col = -1;
  int num_alleles = 0;
  int start_gt_combin_num = -1;
  int end_gt_combin_num = -1;
  char fasta[MAX_FILENAME_LEN];
  AssumptionsOnGraphCleaning assump;
  int ret = parse_genotype_site_argument(arg, list_cols, &num_cols_to_gt, &ref_min_site_col, &num_alleles, &start_gt_combin_num, &end_gt_combin_num, fasta, &assump);
  
  CU_ASSERT(num_cols_to_gt==2);
  CU_ASSERT(list_cols[0]==1);
  CU_ASSERT(list_cols[1]==2);
  CU_ASSERT(ref_min_site_col==3);
  CU_ASSERT(num_alleles==3);
  CU_ASSERT(start_gt_combin_num==2);
  CU_ASSERT(end_gt_combin_num==4);
  CU_ASSERT(strcmp(fasta, "zam/iqbal")==0);
  CU_ASSERT(assump==AssumeUncleaned);






  char arg2[]="*0,1,2,3,4,5[5[5[-1[zam/iqbal[CLEANED"; //means gt all colours except 0,1,2,3,4,5 and that 5 is the ref-minus-site colours, and there are 5 alleles,
                                       //and that we must gt all of the 5 choose 2  = 10 possible genotypes

  int list_cols2[NUMBER_OF_COLOURS];
  int num_cols_to_gt2=0;
  int ref_min_site_col2 = -1;
  int num_alleles2 = 0;
  int start_gt_combin_num2 = -1;
  int end_gt_combin_num2 = -1;
  char fasta2[MAX_FILENAME_LEN];
  AssumptionsOnGraphCleaning assump2;

  
  int ret2 = parse_genotype_site_argument(arg2, list_cols2, &num_cols_to_gt2, &ref_min_site_col2, &num_alleles2, &start_gt_combin_num2, &end_gt_combin_num2, 
					  fasta2, &assump2);
  
  if (NUMBER_OF_COLOURS-6<=0)
    {
      die("We need >6 colours to run the cmdline tests. Recompile.\n");
    }
  CU_ASSERT(num_cols_to_gt2==(NUMBER_OF_COLOURS-6));
  int j;
  for (j=6; j<NUMBER_OF_COLOURS; j++)
    {
      CU_ASSERT(list_cols2[j-6]==j);
    }

  CU_ASSERT(ref_min_site_col2==5);
  CU_ASSERT(num_alleles2==5);
  CU_ASSERT(start_gt_combin_num2==1);
  CU_ASSERT(end_gt_combin_num2==10);
  CU_ASSERT(strcmp(fasta2, "zam/iqbal")==0);
  CU_ASSERT(assump2==AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice);



  char arg3[]="3,4,5[6[3[1,4[zam2/iqbal[CLEANED"; //means gt colours 3,4,5 and that 6 is the ref-minus-site colours, and there are 3 alleles,
                                       //and that we must gt genotypes numbered 1 to 4 out of the 3 choose 2 = 6.

  int list_cols3[NUMBER_OF_COLOURS];
  int num_cols_to_gt3=0;
  int ref_min_site_col3 = -1;
  int num_alleles3 = 0;
  int start_gt_combin_num3 = -1;
  int end_gt_combin_num3 = -1;
  char fasta3[MAX_FILENAME_LEN];
  AssumptionsOnGraphCleaning assump3;

  
  int ret3 = parse_genotype_site_argument(arg3, list_cols3, &num_cols_to_gt3, &ref_min_site_col3, &num_alleles3, &start_gt_combin_num3, &end_gt_combin_num3, fasta3, &assump3);
  
  if (NUMBER_OF_COLOURS-7<=0)
    {
      die("We need >7 colours to run the cmdline tests. Recompile.\n");
    }
  CU_ASSERT(num_cols_to_gt3==3);

  CU_ASSERT(list_cols3[0]==3);
  CU_ASSERT(list_cols3[1]==4);
  CU_ASSERT(list_cols3[2]==5);
  CU_ASSERT(ref_min_site_col3==6);
  CU_ASSERT(num_alleles3==3);
  CU_ASSERT(start_gt_combin_num3==1);
  CU_ASSERT(end_gt_combin_num3==4);
  CU_ASSERT(strcmp(fasta3, "zam2/iqbal")==0);
  CU_ASSERT(assump3==AssumeAnyErrorSeenMustHaveOccurredAtLeastTwice);




}

void test_parse_cmdline_inner_loop_are_basic_variables_correctly_set()
{


  // --COLOUR_LIST
  char* argv1[] = {"cortex", "--colour_list=../data/test/pop_graph/cmd_line/colour_list_1"};
  
  //set up a default cmdline
  CmdLine test1;
  default_opts(&test1);
  char error_string1[LEN_ERROR_STRING];
  int err = parse_cmdline_inner_loop(2, argv1, sizeof(Element), &test1, error_string1);
  CU_ASSERT(err==0);
  CU_ASSERT_STRING_EQUAL(test1.colour_list, "../data/test/pop_graph/cmd_line/colour_list_1");
  CU_ASSERT(test1.input_colours==true);



  //--MULTICOLOUR_BIN

   char* argv2[] = {"cortex", "--multicolour_bin=../data/test/pop_graph/cmd_line/multicolour_bin_1.ctx"};  
  //set up a default cmdline
  CmdLine test2;
  default_opts(&test2);
  char error_string2[LEN_ERROR_STRING];
  int err2 = parse_cmdline_inner_loop(2, argv2, sizeof(Element), &test2, error_string2);
  CU_ASSERT(err2==0);

  CU_ASSERT_STRING_EQUAL(test2.multicolour_bin, "../data/test/pop_graph/cmd_line/multicolour_bin_1.ctx");
  CU_ASSERT(test2.input_multicol_bin==true);


  // SE LIST

  char* argv3[] = {"cortex", "--se_list=../data/test/pop_graph/cmd_line/se_list1"};
  
  //set up a default cmdline
  CmdLine test3;
  default_opts(&test3);

  char error_string3[LEN_ERROR_STRING];
  int err3 = parse_cmdline_inner_loop(2, argv3, sizeof(Element), &test3, error_string3);

  CU_ASSERT_STRING_EQUAL(test3.se_list, "../data/test/pop_graph/cmd_line/se_list1");
  CU_ASSERT(test3.input_seq==true);

  //PE_LIST

  char* argv4[] = {"cortex", "--pe_list=../data/test/pop_graph/cmd_line/pe_list1,../data/test/pop_graph/cmd_line/pe_list2"};
  
  //set up a default cmdline
  CmdLine test4;
  default_opts(&test4);

  char error_string4[LEN_ERROR_STRING];
  int err4 = parse_cmdline_inner_loop(2, argv4, sizeof(Element), &test4, error_string4);

  CU_ASSERT_STRING_EQUAL(test4.pe_list_lh_mates, "../data/test/pop_graph/cmd_line/pe_list1");
  CU_ASSERT_STRING_EQUAL(test4.pe_list_rh_mates, "../data/test/pop_graph/cmd_line/pe_list2");
  CU_ASSERT(test4.input_seq==true);
  

  // KMER SIZE

  char* argv5[] = {"cortex", "--kmer_size=17"};
  
  //set up a default cmdline
  CmdLine test5;
  default_opts(&test5);

  char error_string5[LEN_ERROR_STRING];
  int err5 = parse_cmdline_inner_loop(2, argv5, sizeof(Element), &test5, error_string5);

  CU_ASSERT(test5.kmer_size==17); 


  // HASH TABLE NUMBER OF BUCKETS

  char* argv6[] = {"cortex", "--mem_height=20"};
  
  //set up a default cmdline
  CmdLine test6;
  default_opts(&test6);

  char error_string6[LEN_ERROR_STRING];
  int err6 = parse_cmdline_inner_loop(2, argv6, sizeof(Element), &test6, error_string6);
  
  CU_ASSERT(test6.number_of_buckets_bits==20);


  // DEPTH OF HASH TABLE BUCKETS
  
  char* argv7[] = {"cortex", "--mem_height=100"};
  
  //set up a default cmdline
  CmdLine test7;
  default_opts(&test7);

  char error_string7[LEN_ERROR_STRING];
  int err7 = parse_cmdline_inner_loop(2, argv7, sizeof(Element), &test7, error_string7);
  
  CU_ASSERT(test7.bucket_size==100);
  

  // REF COLOUR

  char* argv8[] = {"cortex", "--ref_colour=0"};
  
  //set up a default cmdline
  CmdLine test8;
  default_opts(&test8);

  char error_string8[LEN_ERROR_STRING];
  int err8 = parse_cmdline_inner_loop(2, argv8, sizeof(Element), &test8, error_string8);
  
  CU_ASSERT(test8.ref_colour==0);
  CU_ASSERT(test8.using_ref==true);

  // REMV PCR DUPS

  char* argv9[] = {"cortex", "--remove_pcr_duplicates"};
  
  //set up a default cmdline
  CmdLine test9;
  default_opts(&test9);

  char error_string9[LEN_ERROR_STRING];
  int err9 = parse_cmdline_inner_loop(2, argv9, sizeof(Element), &test9, error_string9);
  
  CU_ASSERT(test9.remove_pcr_dups==true);


  // CUT HOMOPOLYMERS

  char* argv10[] = {"cortex", "--cut_homopolymers=4"};
  
  //set up a default cmdline
  CmdLine test10;
  default_opts(&test10);

  char error_string10[LEN_ERROR_STRING];
  int err10 = parse_cmdline_inner_loop(2, argv10, sizeof(Element), &test10, error_string10);
  
  CU_ASSERT(test10.cut_homopolymers==true);
  CU_ASSERT(test10.homopolymer_limit==4);

  

  
  // QUALITY SCORE THRESHOLD
  
  char* argv12[] = {"cortex", "--quality_score_threshold=47"};
  
  //set up a default cmdline
  CmdLine test12;
  default_opts(&test12);

  char error_string12[LEN_ERROR_STRING];
  int err12 = parse_cmdline_inner_loop(2, argv12, sizeof(Element), &test12, error_string12);
  
  CU_ASSERT(test12.quality_score_threshold==47);



  // REMOVE SEQ ERRORS USING COVG AND TOPOLOGY OF GRAPH


  char* argv14[] = {"cortex", "--remove_seq_errors"};
  
  //set up a default cmdline
  CmdLine test14;
  default_opts(&test14);

  char error_string14[LEN_ERROR_STRING];
  int err14 = parse_cmdline_inner_loop(2, argv14, sizeof(Element), &test14, error_string14);
  
  CU_ASSERT(test14.remove_seq_errors==true);



  // DUMP BINARY

  char* argv15[] = {"cortex", "--dump_binary=/path/to/nonexistent/file"};
  
  //set up a default cmdline
  CmdLine test15;
  default_opts(&test15);

  char error_string15[LEN_ERROR_STRING];
  int err15 = parse_cmdline_inner_loop(2, argv15, sizeof(Element), &test15, error_string15);
  
  CU_ASSERT(test15.dump_binary==true);
  CU_ASSERT_STRING_EQUAL(test15.output_binary_filename, "/path/to/nonexistent/file");

  


  // DETECT_BUBBLES_1
  
  char* argv17[] = {"cortex", "--detect_bubbles1=1,4,9/2,4"};
  
  //set up a default cmdline
  CmdLine test17;
  default_opts(&test17);

  char error_string17[LEN_ERROR_STRING];
  int err17 = parse_cmdline_inner_loop(2, argv17, sizeof(Element), &test17, error_string17);
  
  CU_ASSERT(test17.detect_bubbles1==true);
  CU_ASSERT(test17.num_colours_in_detect_bubbles1_first_colour_list==3);
  CU_ASSERT(test17.detect_bubbles1_first_colour_list[0]==1);
  CU_ASSERT(test17.detect_bubbles1_first_colour_list[1]==4);
  CU_ASSERT(test17.detect_bubbles1_first_colour_list[2]==9);
  CU_ASSERT(test17.num_colours_in_detect_bubbles1_second_colour_list==2);
  CU_ASSERT(test17.detect_bubbles1_second_colour_list[0]==2 );
  CU_ASSERT(test17.detect_bubbles1_second_colour_list[1]==4 );


  // DETECT_BUBBLES_2
  
  char* argv18[] = {"cortex", "--detect_bubbles2=1,4,9/2,4"};
  
  //set up a default cmdline
  CmdLine test18;
  default_opts(&test18);

  char error_string18[LEN_ERROR_STRING];
  int err18 = parse_cmdline_inner_loop(2, argv18, sizeof(Element), &test18, error_string18);
  
  CU_ASSERT(test18.detect_bubbles2==true);
  CU_ASSERT(test18.num_colours_in_detect_bubbles2_first_colour_list==3);
  CU_ASSERT(test18.detect_bubbles2_first_colour_list[0]==1);
  CU_ASSERT(test18.detect_bubbles2_first_colour_list[1]==4);
  CU_ASSERT(test18.detect_bubbles2_first_colour_list[2]==9);
  CU_ASSERT(test18.num_colours_in_detect_bubbles2_second_colour_list==2);
  CU_ASSERT(test18.detect_bubbles2_second_colour_list[0]==2 );
  CU_ASSERT(test18.detect_bubbles2_second_colour_list[1]==4 );

  
  
  // OUTPUT DETECT_BUBBLES1

  char* argv19[] = {"cortex", "--output_bubbles1=/nonexistent_file_zamzam"};
  
  //set up a default cmdline
  CmdLine test19;
  default_opts(&test19);

  char error_string19[LEN_ERROR_STRING];
  int err19 = parse_cmdline_inner_loop(2, argv19, sizeof(Element), &test19, error_string19);
  
  CU_ASSERT_STRING_EQUAL(test19.output_detect_bubbles1, "/nonexistent_file_zamzam");


  // OUTPUT DETECT_BUBBLES2

  char* argv20[] = {"cortex", "--output_bubbles2=/nonexistent_file_zim"};
  
  //set up a default cmdline
  CmdLine test20;
  default_opts(&test20);

  char error_string20[LEN_ERROR_STRING];
  int err20 = parse_cmdline_inner_loop(2, argv20, sizeof(Element), &test20, error_string20);
  
  CU_ASSERT_STRING_EQUAL(test20.output_detect_bubbles2, "/nonexistent_file_zim");


  // INPUT FORMAT

  char* argv21[] = {"cortex", "--format=CTX"};
  
  //set up a default cmdline
  CmdLine test21;
  default_opts(&test21);

  char error_string21[LEN_ERROR_STRING];
  int err21 = parse_cmdline_inner_loop(2, argv21, sizeof(Element), &test21, error_string21);
  
  CU_ASSERT(test21.format_of_input_seq==CTX);


  //MAX READ LENGTH

  char* argv22[] = {"cortex", "--max_read_len=1047"};
  
  //set up a default cmdline
  CmdLine test22;
  default_opts(&test22);

  char error_string22[LEN_ERROR_STRING];
  int err22 = parse_cmdline_inner_loop(2, argv22, sizeof(Element), &test22, error_string22);
  
  CU_ASSERT(test22.max_read_length==1047);



  // PRINT COLOUR COVERAGES

  char* argv23[] = {"cortex", "--print_colour_coverages"};
  
  //set up a default cmdline
  CmdLine test23;
  default_opts(&test23);

  char error_string23[LEN_ERROR_STRING];
  int err23 = parse_cmdline_inner_loop(2, argv23, sizeof(Element), &test23, error_string23);

  CU_ASSERT(err23==0);
  CU_ASSERT(test23.print_colour_coverages==true);
  
}




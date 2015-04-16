/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 *
 */
/*
  run_cmd_line_tests.c
*/

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "cmd_line.h"
#include "test_cmd_line.h"
//#include "test_dB_graph.h"

int  main()
{

  CU_pSuite pCortexVar_CmdLine_Test_Suite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();


  /* add a suite to the registry */
  pCortexVar_CmdLine_Test_Suite = CU_add_suite("Tests of Command-Line of cortex_var - 1", NULL, NULL);
  if (NULL == pCortexVar_CmdLine_Test_Suite) {
    CU_cleanup_registry();
    return CU_get_error();
  }


  /* add the tests to the suites */

  if (NULL == CU_add_test(pCortexVar_CmdLine_Test_Suite, "Test utility function for getting colours from input format", test_get_numbers_from_comma_sep_list)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  

  if (NULL == CU_add_test(pCortexVar_CmdLine_Test_Suite, "Test utility function for parsing args (lists of colours separated by /) for detect_bubbles", test_parse_colourinfo_argument)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pCortexVar_CmdLine_Test_Suite, "Test utility function for parsing args (lists of colours separated by ,) for path divergence call specification", test_parse_commasep_list)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pCortexVar_CmdLine_Test_Suite, "Test utility function for parsing  argument of --genotype_site", test_parse_genotype_site_argument)) {
    CU_cleanup_registry();
    return CU_get_error();
  }



  if (NULL == CU_add_test(pCortexVar_CmdLine_Test_Suite, "Test that when given good input, internal variables are correctly set", test_parse_cmdline_inner_loop_are_basic_variables_correctly_set )) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  //  if (NULL == CU_add_test(pCortexVar_CmdLine_Test_Suite, "Test that when given good input, internal variables are correctly set 2", test_parse_cmdline_inner_loop_are_basic_variables_correctly_set_2 )) {
  //  CU_cleanup_registry();
  //  return CU_get_error();
  //}


  printf("Start 6\n");



 

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}


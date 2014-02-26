
//#include "test_dB_graph.h"
#include "test_build.h"
#include "test_genotyping_known.h"
#include "test_gene_presence.h"
#include <CUnit.h>
#include <Basic.h>

int  main()
{

  CU_pSuite pPopGraphSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();

  /* add a suite to the registry */
  pPopGraphSuite = CU_add_suite("Test main functionality of myKrobe predictor", NULL, NULL);
  if (NULL == pPopGraphSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suites */



  /*  if (NULL == CU_add_test(pPopGraphSuite, "Test building unclean graph", test_build_unclean_graph)) {
    CU_cleanup_registry();
    return CU_get_error();
    }*/

  if (NULL == CU_add_test(pPopGraphSuite, "Test getting coverage info on resistance/susceptibility alleles", test_get_next_mutation_allele_info)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, "Test getting coverage info on a gene (for gene presence testing)", test_get_next_gene_info)) {
    CU_cleanup_registry();
    return CU_get_error();
  }


  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}


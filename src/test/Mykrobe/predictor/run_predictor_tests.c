
//#include "test_dB_graph.h"
#include "test_build.h"
#include "test_genotyping_known.h"
#include "test_gene_presence.h"
#include "test_species_prediction.h"
#include "test_gene_presence_models.h"
#include "test_mut_models.h"
#include <CUnit.h>
#include <Basic.h>

int  main()
{

  CU_pSuite pPopGraphSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();

  /* add a suite to the registry */
  pPopGraphSuite = CU_add_suite("Test main functionality of Mykrobe predictor", NULL, NULL);
  if (NULL == pPopGraphSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suites */



  /*  if (NULL == CU_add_test(pPopGraphSuite, "Test building unclean graph", test_build_unclean_graph)) {
    CU_cleanup_registry();
    return CU_get_error();
    }*/

  // if (NULL == CU_add_test(pPopGraphSuite, 
		//  	  "Test getting coverage info on resistance/susceptibility alleles", 
		//  	  test_get_next_var_on_background)) 
  //   {
  //     CU_cleanup_registry();
  //     return CU_get_error();
  //   }

  //  if (NULL == CU_add_test(pPopGraphSuite, 
		//  	  "Test log likelihoods/models for clonal susceptible infection", 
		//  	  test_mutation_model_log_likelihoods_1)) {
  //    CU_cleanup_registry();
  //    return CU_get_error();
  //  }

  //  if (NULL == CU_add_test(pPopGraphSuite, 
		//  	  "Test log likelihoods/models for clonal resistant infection", 
		//  	  test_mutation_model_log_likelihoods_2)) {
  //    CU_cleanup_registry();
  //    return CU_get_error();
  //  }
  if (NULL == CU_add_test(pPopGraphSuite, 
			  "Test gene presence models ", 
			  test_resistotype_gene)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test gene presence minor models ", 
       test_resistotype_minor_gene)) {
    CU_cleanup_registry();
    return CU_get_error();
  }  
  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test gene presence unsure models ", 
       test_resistotype_unsure_gene)) {
    CU_cleanup_registry();
    return CU_get_error();
  }  

  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test gene presence unsure models 2 ", 
       test_resistotype_unsure_gene_2)) {
    CU_cleanup_registry();
    return CU_get_error();
  }  


  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test gene presence high CN models ", 
       test_resistotype_gene_at_high_CN)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test gene not present models ", 
       test_resistotype_gene_S)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test gene low coverage ", 
       test_low_coverage_genes)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test gene low coverage ont", 
       test_low_coverage_ont_genes)) {
    CU_cleanup_registry();
    return CU_get_error();
  }


 // Mutation Models
  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test mutation not present models ", 
       test_mutation_S)) {
    CU_cleanup_registry();
    return CU_get_error();
  }
  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test mutation R models ", 
       test_mutation_R)) {
    CU_cleanup_registry();
    return CU_get_error();
  }  
  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test mutation minor models ", 
       test_mutation_r)) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pPopGraphSuite, 
       "Test mutation custom case ", 
       test_mutation_custom_1)) {
    CU_cleanup_registry();
    return CU_get_error();
  }



  // if (NULL == CU_add_test(pPopGraphSuite, 
		// 	  "Test getting coverage info on a gene (for gene presence testing)", 
		// 	  test_get_next_gene_info)) {
  //   CU_cleanup_registry();
  //   return CU_get_error();
  // }
  // run to get comments before test
  // test_get_species_info();
  // if (NULL == CU_add_test(pPopGraphSuite, "Test assigning a species to a sample", test_get_species_info)) {
  //   CU_cleanup_registry();
  //   return CU_get_error();
  // }

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}


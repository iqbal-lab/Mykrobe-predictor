/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 *
  run_hash_table_tests.c
*/

#include <CUnit.h>
#include <Basic.h>

// cortex_var headers
#include "test_hash.h"

int  main()
{

  CU_pSuite pSuite = NULL;

  /* initialize the CUnit test registry */
  if (CUE_SUCCESS!=CU_initialize_registry())
    return CU_get_error();
  
  /* add a suite to the registry */
  pSuite = CU_add_suite("Suite_1", NULL, NULL);
  if (NULL == pSuite) {
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* add the tests to the suite */

  
  if (NULL == CU_add_test(pSuite, "test hash_table_find_or_insert using Element from graph/",  test_hash_table_find_or_insert)){
    CU_cleanup_registry();
    return CU_get_error();
  }

  if (NULL == CU_add_test(pSuite, "test hash_table_apply_or_insert using Element from graph/. This test takes a little while - do not panic if it takes a minute or two",  test_hash_table_apply_or_insert)){
    CU_cleanup_registry();
    return CU_get_error();
  }

  /* Run all tests using the CUnit Basic interface */
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
 
  CU_cleanup_registry();
  return CU_get_error();


}


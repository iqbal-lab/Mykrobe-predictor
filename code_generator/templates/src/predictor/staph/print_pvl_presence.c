void print_{{gene|lower}}_presence(dBGraph* db_graph,
			int (*file_reader)(FILE * fp, 
					   Sequence * seq, 
					   int max_read_length, 
					   boolean new_entry, 
					   boolean * full_entry),
			ReadingUtils* rutils,
			GeneInfo* tmp_gi,
			Troolean (*func)(dBGraph* db_graph,
					int (*file_reader)(FILE * fp, 
							   Sequence * seq, 
							   int max_read_length, 
							   boolean new_entry, 
							   boolean * full_entry),
					ReadingUtils* rutils,
					GeneInfo* tmp_gi,
					StrBuf* install_dir),
			StrBuf* install_dir, OutputFormat format)
{

  Troolean result = is_{{gene|lower}}_positive(db_graph, file_reader, rutils, tmp_gi, install_dir);
  
  if (format==Stdout)
    {
      printf("{{gene|lower}}\t");
      if (result==true)
	{
	  printf("positive\n");
	}
      else
	{
	  printf("negative\n");
	}
    }
  else
    {
      print_json_virulence_start();
      if (result==true)
	{
	  print_json_item("{{gene|lower}}", "positive", {% if loop.last %}true {% else false%});
	}
      else
	{
	  print_json_item("{{gene|lower}}", "negative", {% if loop.last %}true {% else false%});
	}
      print_json_virulence_end();
    }
}
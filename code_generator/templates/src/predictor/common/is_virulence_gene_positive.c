Troolean is_{{gene| lower}}_positive(dBGraph* db_graph,
			   int (*file_reader)(FILE * fp, 
					      Sequence * seq, 
					      int max_read_length, 
					      boolean new_entry, 
					      boolean * full_entry),
			ReadingUtils* rutils,
			GeneInfo* tmp_gi,
			StrBuf* install_dir)
			   

{
  StrBuf* fa = strbuf_create(install_dir->buff);
  strbuf_append_str(fa, "data/{{selfer.species}}/virulence/{{gene}}.fa");

  FILE* fp = fopen(fa->buff, "r");
  if (fp==NULL)
    {
      die("Cannot open %s\n", fa->buff);
    }
  int num=1;
  boolean is_pos=false;
  while (num>0)
    {
      num = get_next_gene_info(fp,
			       db_graph,
			       tmp_gi,
			       rutils->seq,
			       rutils->kmer_window,
			       file_reader,
			       rutils->array_nodes,
			       rutils->array_or,
			       rutils->working_ca,
			       MAX_LEN_GENE);
      if (tmp_gi->percent_nonzero > MIN_PERC_COVG_VIRULENCE)
	{
	  is_pos=true;
	}
    }
  fclose(fp);
  strbuf_free(fa);
  return is_pos;

}
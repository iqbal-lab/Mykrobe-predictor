GeneMutationGene map_gene_name_str_to_genename(StrBuf* name)
{
  {% for gene in selfer.all_mut_genes %}
  {% if loop.first %}if{% else %}else if{% endif %}(strcmp(name->buff, "{{gene}}")==0)
    {
      return {{gene}};
    }
  {% endfor %}
  else
    {
      die("Failed to parse gene name - got %s\n", name->buff);
    }

}


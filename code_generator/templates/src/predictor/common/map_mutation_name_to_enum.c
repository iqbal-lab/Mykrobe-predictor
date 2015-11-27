
KnownMutation map_mutation_name_to_enum(StrBuf* sbuf, GeneMutationGene gene)
{
  {% for mut in selfer.all_mutations %}
  
  {% if loop.first %}if{% else %}else if{% endif %} ( (strcmp(sbuf->buff, "{{mut.split('_')[-1]}}")==0) && (gene=={{mut.split('_')[0]}}) )
    {
      return {{mut}};
    } 
  {% endfor %}
  else 
    {
      die("Parsing error - unknown mutation %s for gene %d\n", sbuf->buff, (int) gene);
      return NotSpecified;
    } 
}
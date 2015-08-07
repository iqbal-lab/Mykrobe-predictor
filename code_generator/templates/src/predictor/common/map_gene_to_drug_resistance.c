char* map_gene_to_drug_resistance(GenePresenceGene gene)
{
{% if selfer.genes %}
   switch (gene) 
   {
    {% for gene_enum in selfer.genes %}
    case {{gene_enum}} : return "{{",".join(selfer.gene_enum_to_drug_name[gene_enum])}}";
    {% endfor %}
   }
{% endif %}	   
   return "unknown";
}

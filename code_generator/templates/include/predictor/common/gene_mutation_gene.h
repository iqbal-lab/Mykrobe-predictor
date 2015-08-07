  typedef enum
  {
    {% for gene in selfer.all_mut_genes %}
    {{gene}}={{loop.index0}},
    {% endfor %}
    Unknown = {{selfer.all_mut_genes|length}}
  }GeneMutationGene;
  #define NUM_KNOWN_GENES {{selfer.all_mut_genes|length + 1 }}
/*
 * Copyright 2015 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
  species.c
*/
// system headers
#include <stdlib.h>
#include <limits.h>

// third party headers
#include <string_buffer.h>

#include "build.h" 
#include "maths.h" 
#include "element.h"
#include "seq.h"
#include "open_hash/hash_table.h"
#include "dB_graph.h"
#include "species.h"
#include "gene_presence.h"
#include "genotyping_known.h"

{% for phylogroup in selfer.phylo_groups %}
  {% include 'src/predictor/common/phylo/map_phylo_group_enum_to_str.c' %}
  {% include 'src/predictor/common/phylo/load_all_phylo_group_file_paths.c' %}
  {% include 'src/predictor/common/phylo/phylo_group_threshold.c' %}
  {% if selfer.species == "staph" and phylogroup.name in ["species","phylo_group"] %}
    {% else %}
    {% include 'src/predictor/common/phylo/print_json_phylo_group.c' %}
  {% endif %}
  {% include 'src/predictor/common/phylo/get_ith_phylo_group_name.c' %}
{% endfor %}


{% if not selfer.species == "staph" %}
SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last)

{
{% for phylogroup in selfer.phylo_groups %}
  char* get_ith_{{phylogroup.name}}_name(CovgInfo* covg_info, int i);
  StrBuf* {{phylogroup.name}}_file_paths[NUM_{{phylogroup.enum}}];
  load_all_{{phylogroup.name}}_file_paths({{phylogroup.name}}_file_paths,install_dir);
  CovgInfo* {{phylogroup.name}}_covg_info = get_coverage_info(db_graph,
                                                  {{phylogroup.name}}_file_paths,
                                                  max_branch_len,NUM_{{phylogroup.enum}},
                                                  ignore_first,ignore_last,
                                                  {{phylogroup.name}}_threshold);  

{% endfor %}


  SpeciesInfo* species_info=(SpeciesInfo *)malloc(sizeof(SpeciesInfo)); 
  {% for phylogroup in selfer.phylo_groups %}  
  species_info->{{phylogroup.name}}_covg_info = {{phylogroup.name}}_covg_info;
  {% endfor %}

  // update_phylo_group_presence_and_coverage_from_species(species_info);

  return species_info;
}
{% endif %}


{% for phylogroup in selfer.phylo_groups %}

void print_json_{{phylogroup.name}}_start()
{
  printf("\t\t\"{{phylogroup.name}}\": {\n");
}
{% endfor %}

void print_json_phylogenetics(SpeciesInfo* species_info){
    print_json_phylogenetics_start();

  {% for phylogroup in selfer.phylo_groups %}
    print_json_{{phylogroup.name}}(species_info);
  {% endfor %}
  {% if selfer.species == "staph" %}
      print_json_lineage(species_info);
  {% endif %}

    print_json_phylogenetics_end();  
}

{% block extra %}
{% endblock %}

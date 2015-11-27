 typedef enum
  {
    {% for mut in selfer.all_mutations %}
    {{mut}}={{loop.index0}},
    {% endfor %}
    NotSpecified = {{selfer.all_mutations | length}}
  } KnownMutation;
#define NUM_KNOWN_MUTATIONS {{selfer.all_mutations | length + 1}}


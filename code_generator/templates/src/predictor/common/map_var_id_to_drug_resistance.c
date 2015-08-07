{% if selfer.all_mutations %}
char* map_var_id_to_drug_resistance(KnownMutation km)
{
   switch (km) 
   {
    {% for drug in selfer.drugs %}
        {% for mut in drug.mut_list %}
        case {{mut}}  : return "{{drug}}";
        {% endfor %}
    {% endfor %}
     case  NotSpecified : return "NotSpecified";      
   }
   return  "NotSpecified";
}
{% else %}
char* map_var_id_to_drug_resistance(KnownMutation km)
{
   return  "NotSpecified";
}

{% endif %}
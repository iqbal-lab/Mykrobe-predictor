{% if selfer.all_mutations %}
	const char* map_enum_to_mutation_name(KnownMutation km)
	{
	   switch (km) 
	   {
	    {% for mut in selfer.all_mutations %}
	     case {{mut}}  : return "{{mut}}";
	     {% endfor %}
	     case  NotSpecified : return "NotSpecified";      
	   }
	   return  "NotSpecified";
	}
{% else %}
	const char* map_enum_to_mutation_name(KnownMutation km)
	{
	   return  "NotSpecified";
	}
{% endif %}
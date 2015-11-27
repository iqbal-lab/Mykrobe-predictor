typedef enum 
 {
   NoDrug=0,
   {% for drug in selfer.drugs %}
   {{drug.name}}={{loop.index}}{% if not loop.last %},{% endif %}
   {% endfor %}
  } Antibiotic;
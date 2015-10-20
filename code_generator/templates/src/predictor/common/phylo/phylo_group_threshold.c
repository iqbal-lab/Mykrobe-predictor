void {{phylogroup.name}}_threshold(int* thresholds){
	{% for taxon in phylogroup.taxons %}
	  thresholds[{{taxon.enum}}] = {{taxon.threshold}};
	{% endfor %}
}
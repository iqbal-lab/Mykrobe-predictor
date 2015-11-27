void load_all_{{phylogroup.name}}_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
	{% for taxon in phylogroup.taxons %}
  panel_file_paths[{{taxon.enum}}] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[{{taxon.enum}}], "data/{{selfer.species}}/phylo/{{phylogroup.name}}/{{taxon.filename}}" );
  {% endfor %}
}
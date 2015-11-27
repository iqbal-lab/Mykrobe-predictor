char* get_ith_{{phylogroup.name}}_name(CovgInfo* covg_info, int i)
{
  {{phylogroup.enum}} {{phylogroup.name}};
  StrBuf* {{phylogroup.name}}_name = strbuf_new(); 
  {{phylogroup.name}} = get_ith_present_panel( covg_info, i);
  map_{{phylogroup.name}}_enum_to_str({{phylogroup.name}}, {{phylogroup.name}}_name);
  return {{phylogroup.name}}_name->buff;
}
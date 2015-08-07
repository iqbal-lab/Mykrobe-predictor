void map_antibiotic_enum_to_str(Antibiotic ab, StrBuf* name)
{
  if (ab==NoDrug)
    {
      strbuf_reset(name);
      strbuf_append_str(name, "NoDrug");
    }
  {% for drug_enum in selfer.drugs %}
  else if (ab=={{drug_enum}})
    {
      strbuf_reset(name);
      strbuf_append_str(name, "{{drug_enum}}");
    }
  {% endfor %}
  else
    {
      die("Impossible - compiler should not allow this\n");
    }
}
void map_{{phylogroup.name}}_enum_to_str({{ phylogroup.enum }} sp, StrBuf* sbuf)
{
  {% for taxon in phylogroup.taxons %}
    {% if loop.first %} if{% else %} else if{% endif %}(sp=={{taxon.enum}}){
        strbuf_reset(sbuf);
        strbuf_append_str(sbuf, "{{taxon.name}}");
      }
  {% endfor %}
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value %i \n",sp );
    }
}
/*
 * Copyright 2014 Zamin Iqbal (zam@well.ox.ac.uk)
 * 
 *
 * **********************************************************************
 *
 * This file is part of Mykrobe.
 *
 * Mykrobe is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Mykrobe is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Mykrobe.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
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

void map_phylo_group_enum_to_str(PhyloGroup sp, StrBuf* sbuf)
{
  if(sp==CoagPos){
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Staphylococcus aureus");
  }
  else if(sp==CoagNeg){
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Coagulase-Negative Staphylococcus");
  }  
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value %i \n",sp );
    }
}

void map_species_enum_to_str(Species species, StrBuf* sbuf)
{
  if (species == Achromobacter_piechaudii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Achromobacter_piechaudii");
  }
else if (species == Achromobacter_xylosoxidans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Achromobacter_xylosoxidans");
  }
else if (species == Acidaminococcus_fermentans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acidaminococcus_fermentans");
  }
else if (species == Acidithiobacillus_caldus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acidithiobacillus_caldus");
  }
else if (species == Acidithiobacillus_ferrooxidans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acidithiobacillus_ferrooxidans");
  }
else if (species == Acidovorax_avenae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acidovorax_avenae");
  }
else if (species == Acidovorax_delafieldii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acidovorax_delafieldii");
  }
else if (species == Acinetobacter_baumannii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acinetobacter_baumannii");
  }
else if (species == Acinetobacter_calcoaceticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acinetobacter_calcoaceticus");
  }
else if (species == Acinetobacter_haemolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acinetobacter_haemolyticus");
  }
else if (species == Acinetobacter_johnsonii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acinetobacter_johnsonii");
  }
else if (species == Acinetobacter_junii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acinetobacter_junii");
  }
else if (species == Acinetobacter_lwoffii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acinetobacter_lwoffii");
  }
else if (species == Acinetobacter_radioresistens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Acinetobacter_radioresistens");
  }
else if (species == Actinobacillus_minor)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Actinobacillus_minor");
  }
else if (species == Actinobacillus_pleuropneumoniae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Actinobacillus_pleuropneumoniae");
  }
else if (species == Actinobacillus_succinogenes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Actinobacillus_succinogenes");
  }
else if (species == Actinobacillus_ureae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Actinobacillus_ureae");
  }
else if (species == Actinomyces_coleocanis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Actinomyces_coleocanis");
  }
else if (species == Actinomyces_odontolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Actinomyces_odontolyticus");
  }
else if (species == Actinomyces_oris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Actinomyces_oris");
  }
else if (species == Actinomyces_urogenitalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Actinomyces_urogenitalis");
  }
else if (species == Actinomyces_viscosus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Actinomyces_viscosus");
  }
else if (species == Aeromonas_hydrophila)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Aeromonas_hydrophila");
  }
else if (species == Aeromonas_salmonicida)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Aeromonas_salmonicida");
  }
else if (species == Aggregatibacter_actinomycetemcomitans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Aggregatibacter_actinomycetemcomitans");
  }
else if (species == Aggregatibacter_aphrophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Aggregatibacter_aphrophilus");
  }
else if (species == Aggregatibacter_segnis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Aggregatibacter_segnis");
  }
else if (species == Agrobacterium_tumefaciens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Agrobacterium_tumefaciens");
  }
else if (species == Agrobacterium_vitis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Agrobacterium_vitis");
  }
else if (species == Alcanivorax_borkumensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Alcanivorax_borkumensis");
  }
else if (species == Alicyclobacillus_acidocaldarius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Alicyclobacillus_acidocaldarius");
  }
else if (species == Alicyclobacillus_Alicyclobacillus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Alicyclobacillus_Alicyclobacillus");
  }
else if (species == Aliivibrio_fischeri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Aliivibrio_fischeri");
  }
else if (species == Aliivibrio_salmonicida)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Aliivibrio_salmonicida");
  }
else if (species == Alistipes_putredinis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Alistipes_putredinis");
  }
else if (species == Alistipes_shahii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Alistipes_shahii");
  }
else if (species == Alkaliphilus_metalliredigens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Alkaliphilus_metalliredigens");
  }
else if (species == Alkaliphilus_oremlandii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Alkaliphilus_oremlandii");
  }
else if (species == Anabaena_azollae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anabaena_azollae'");
  }
else if (species == Anabaena_variabilis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anabaena_variabilis");
  }
else if (species == Anaerococcus_hydrogenalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaerococcus_hydrogenalis");
  }
else if (species == Anaerococcus_lactolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaerococcus_lactolyticus");
  }
else if (species == Anaerococcus_prevotii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaerococcus_prevotii");
  }
else if (species == Anaerococcus_tetradius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaerococcus_tetradius");
  }
else if (species == Anaerococcus_vaginalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaerococcus_vaginalis");
  }
else if (species == Anaeromyxobacter_dehalogenans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaeromyxobacter_dehalogenans");
  }
else if (species == Anaerostipes_caccae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaerostipes_caccae");
  }
else if (species == Anaplasma_centrale)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaplasma_centrale");
  }
else if (species == Anaplasma_marginale)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaplasma_marginale");
  }
else if (species == Anaplasma_phagocytophilum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Anaplasma_phagocytophilum");
  }
else if (species == Archaeoglobus_fulgidus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Archaeoglobus_fulgidus");
  }
else if (species == Archaeoglobus_profundus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Archaeoglobus_profundus");
  }
else if (species == Arcobacter_butzleri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Arcobacter_butzleri");
  }
else if (species == Arcobacter_nitrofigilis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Arcobacter_nitrofigilis");
  }
else if (species == Arthrobacter_arilaitensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Arthrobacter_arilaitensis");
  }
else if (species == Arthrobacter_aurescens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Arthrobacter_aurescens");
  }
else if (species == Arthrobacter_chlorophenolicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Arthrobacter_chlorophenolicus");
  }
else if (species == Arthrobacter_phenanthrenivorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Arthrobacter_phenanthrenivorans");
  }
else if (species == Arthrospira_maxima)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Arthrospira_maxima");
  }
else if (species == Arthrospira_platensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Arthrospira_platensis");
  }
else if (species == Atopobium_parvulum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Atopobium_parvulum");
  }
else if (species == Atopobium_rimae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Atopobium_rimae");
  }
else if (species == Atopobium_vaginae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Atopobium_vaginae");
  }
else if (species == Bacillus_amyloliquefaciens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_amyloliquefaciens");
  }
else if (species == Bacillus_anthracis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_anthracis");
  }
else if (species == Bacillus_atrophaeus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_atrophaeus");
  }
else if (species == Bacillus_cellulosilyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_cellulosilyticus");
  }
else if (species == Bacillus_cereus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_cereus");
  }
else if (species == Bacillus_clausii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_clausii");
  }
else if (species == Bacillus_coagulans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_coagulans");
  }
else if (species == Bacillus_coahuilensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_coahuilensis");
  }
else if (species == Bacillus_halodurans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_halodurans");
  }
else if (species == Bacillus_licheniformis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_licheniformis");
  }
else if (species == Bacillus_megaterium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_megaterium");
  }
else if (species == Bacillus_mycoides)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_mycoides");
  }
else if (species == Bacillus_pseudofirmus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_pseudofirmus");
  }
else if (species == Bacillus_pseudomycoides)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_pseudomycoides");
  }
else if (species == Bacillus_pumilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_pumilus");
  }
else if (species == Bacillus_selenitireducens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_selenitireducens");
  }
else if (species == Bacillus_subtilis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_subtilis");
  }
else if (species == Bacillus_thuringiensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_thuringiensis");
  }
else if (species == Bacillus_weihenstephanensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacillus_weihenstephanensis");
  }
else if (species == Bacteroides_caccae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_caccae");
  }
else if (species == Bacteroides_cellulosilyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_cellulosilyticus");
  }
else if (species == Bacteroides_coprocola)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_coprocola");
  }
else if (species == Bacteroides_coprophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_coprophilus");
  }
else if (species == Bacteroides_dorei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_dorei");
  }
else if (species == Bacteroides_eggerthii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_eggerthii");
  }
else if (species == Bacteroides_finegoldii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_finegoldii");
  }
else if (species == Bacteroides_fragilis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_fragilis");
  }
else if (species == Bacteroides_helcogenes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_helcogenes");
  }
else if (species == Bacteroides_intestinalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_intestinalis");
  }
else if (species == Bacteroides_ovatus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_ovatus");
  }
else if (species == Bacteroides_pectinophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_pectinophilus");
  }
else if (species == Bacteroides_plebeius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_plebeius");
  }
else if (species == Bacteroides_salanitronis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_salanitronis");
  }
else if (species == Bacteroides_stercoris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_stercoris");
  }
else if (species == Bacteroides_thetaiotaomicron)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_thetaiotaomicron");
  }
else if (species == Bacteroides_uniformis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_uniformis");
  }
else if (species == Bacteroides_vulgatus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_vulgatus");
  }
else if (species == Bacteroides_xylanisolvens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bacteroides_xylanisolvens");
  }
else if (species == Bartonella_bacilliformis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bartonella_bacilliformis");
  }
else if (species == Bartonella_clarridgeiae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bartonella_clarridgeiae");
  }
else if (species == Bartonella_grahamii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bartonella_grahamii");
  }
else if (species == Bartonella_henselae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bartonella_henselae");
  }
else if (species == Bartonella_quintana)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bartonella_quintana");
  }
else if (species == Bartonella_tribocorum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bartonella_tribocorum");
  }
else if (species == Bifidobacterium_adolescentis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_adolescentis");
  }
else if (species == Bifidobacterium_angulatum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_angulatum");
  }
else if (species == Bifidobacterium_animalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_animalis");
  }
else if (species == Bifidobacterium_bifidum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_bifidum");
  }
else if (species == Bifidobacterium_breve)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_breve");
  }
else if (species == Bifidobacterium_catenulatum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_catenulatum");
  }
else if (species == Bifidobacterium_dentium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_dentium");
  }
else if (species == Bifidobacterium_gallicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_gallicum");
  }
else if (species == Bifidobacterium_longum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_longum");
  }
else if (species == Bifidobacterium_pseudocatenulatum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bifidobacterium_pseudocatenulatum");
  }
else if (species == Blautia_hansenii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Blautia_hansenii");
  }
else if (species == Blautia_hydrogenotrophica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Blautia_hydrogenotrophica");
  }
else if (species == Bordetella_avium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bordetella_avium");
  }
else if (species == Bordetella_bronchiseptica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bordetella_bronchiseptica");
  }
else if (species == Bordetella_parapertussis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bordetella_parapertussis");
  }
else if (species == Bordetella_pertussis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bordetella_pertussis");
  }
else if (species == Bordetella_petrii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bordetella_petrii");
  }
else if (species == Borrelia_afzelii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Borrelia_afzelii");
  }
else if (species == Borrelia_burgdorferi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Borrelia_burgdorferi");
  }
else if (species == Borrelia_duttonii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Borrelia_duttonii");
  }
else if (species == Borrelia_garinii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Borrelia_garinii");
  }
else if (species == Borrelia_hermsii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Borrelia_hermsii");
  }
else if (species == Borrelia_recurrentis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Borrelia_recurrentis");
  }
else if (species == Borrelia_spielmanii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Borrelia_spielmanii");
  }
else if (species == Borrelia_turicatae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Borrelia_turicatae");
  }
else if (species == Borrelia_valaisiana)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Borrelia_valaisiana");
  }
else if (species == Brachyspira_hyodysenteriae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brachyspira_hyodysenteriae");
  }
else if (species == Brachyspira_murdochii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brachyspira_murdochii");
  }
else if (species == Brachyspira_pilosicoli)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brachyspira_pilosicoli");
  }
else if (species == Bradyrhizobium_japonicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Bradyrhizobium_japonicum");
  }
else if (species == Brevibacterium_linens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brevibacterium_linens");
  }
else if (species == Brevibacterium_mcbrellneri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brevibacterium_mcbrellneri");
  }
else if (species == Brevundimonas_subvibrioides)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brevundimonas_subvibrioides");
  }
else if (species == Brucella_abortus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brucella_abortus");
  }
else if (species == Brucella_canis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brucella_canis");
  }
else if (species == Brucella_ceti)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brucella_ceti");
  }
else if (species == Brucella_melitensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brucella_melitensis");
  }
else if (species == Brucella_microti)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brucella_microti");
  }
else if (species == Brucella_neotomae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brucella_neotomae");
  }
else if (species == Brucella_ovis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brucella_ovis");
  }
else if (species == Brucella_pinnipedialis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brucella_pinnipedialis");
  }
else if (species == Brucella_suis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Brucella_suis");
  }
else if (species == Burkholderia_ambifaria)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_ambifaria");
  }
else if (species == Burkholderia_cenocepacia)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_cenocepacia");
  }
else if (species == Burkholderia_cepacia)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_cepacia");
  }
else if (species == Burkholderia_dolosa)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_dolosa");
  }
else if (species == Burkholderia_glumae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_glumae");
  }
else if (species == Burkholderia_graminis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_graminis");
  }
else if (species == Burkholderia_mallei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_mallei");
  }
else if (species == Burkholderia_multivorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_multivorans");
  }
else if (species == Burkholderia_oklahomensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_oklahomensis");
  }
else if (species == Burkholderia_phymatum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_phymatum");
  }
else if (species == Burkholderia_phytofirmans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_phytofirmans");
  }
else if (species == Burkholderia_pseudomallei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_pseudomallei");
  }
else if (species == Burkholderia_rhizoxinica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_rhizoxinica");
  }
else if (species == Burkholderia_thailandensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_thailandensis");
  }
else if (species == Burkholderia_ubonensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_ubonensis");
  }
else if (species == Burkholderia_vietnamiensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_vietnamiensis");
  }
else if (species == Burkholderia_xenovorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Burkholderia_xenovorans");
  }
else if (species == Butyrivibrio_crossotus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Butyrivibrio_crossotus");
  }
else if (species == Butyrivibrio_fibrisolvens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Butyrivibrio_fibrisolvens");
  }
else if (species == Butyrivibrio_proteoclasticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Butyrivibrio_proteoclasticus");
  }
else if (species == Caldanaerobacter_pacificum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldanaerobacter_pacificum");
  }
else if (species == Caldanaerobacter_tengcongensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldanaerobacter_tengcongensis");
  }
else if (species == Caldicellulosiruptor_bescii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldicellulosiruptor_bescii");
  }
else if (species == Caldicellulosiruptor_hydrothermalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldicellulosiruptor_hydrothermalis");
  }
else if (species == Caldicellulosiruptor_kristjanssonii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldicellulosiruptor_kristjanssonii");
  }
else if (species == Caldicellulosiruptor_kronotskyensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldicellulosiruptor_kronotskyensis");
  }
else if (species == Caldicellulosiruptor_lactoaceticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldicellulosiruptor_lactoaceticus");
  }
else if (species == Caldicellulosiruptor_obsidiansis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldicellulosiruptor_obsidiansis");
  }
else if (species == Caldicellulosiruptor_owensensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldicellulosiruptor_owensensis");
  }
else if (species == Caldicellulosiruptor_saccharolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caldicellulosiruptor_saccharolyticus");
  }
else if (species == Campylobacter_coli)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_coli");
  }
else if (species == Campylobacter_concisus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_concisus");
  }
else if (species == Campylobacter_curvus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_curvus");
  }
else if (species == Campylobacter_fetus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_fetus");
  }
else if (species == Campylobacter_gracilis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_gracilis");
  }
else if (species == Campylobacter_hominis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_hominis");
  }
else if (species == Campylobacter_jejuni)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_jejuni");
  }
else if (species == Campylobacter_lari)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_lari");
  }
else if (species == Campylobacter_rectus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_rectus");
  }
else if (species == Campylobacter_showae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_showae");
  }
else if (species == Campylobacter_upsaliensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Campylobacter_upsaliensis");
  }
else if (species == Candidatus_Blochmannia_floridanus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Candidatus_Blochmannia_floridanus");
  }
else if (species == Candidatus_Blochmannia_pennsylvanicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Candidatus_Blochmannia_pennsylvanicus");
  }
else if (species == Candidatus_Pelagibacter_ubique)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Candidatus_Pelagibacter_ubique");
  }
else if (species == Candidatus_Phytoplasma_yellows)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Candidatus_Phytoplasma_yellows");
  }
else if (species == Candidatus_Sulcia_muelleri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Candidatus_Sulcia_muelleri");
  }
else if (species == Capnocytophaga_gingivalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Capnocytophaga_gingivalis");
  }
else if (species == Capnocytophaga_ochracea)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Capnocytophaga_ochracea");
  }
else if (species == Capnocytophaga_sputigena)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Capnocytophaga_sputigena");
  }
else if (species == Caulobacter_crescentus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caulobacter_crescentus");
  }
else if (species == Caulobacter_segnis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Caulobacter_segnis");
  }
else if (species == Cellulophaga_algicola)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Cellulophaga_algicola");
  }
else if (species == Cellulophaga_lytica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Cellulophaga_lytica");
  }
else if (species == Chlamydia_muridarum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlamydia_muridarum");
  }
else if (species == Chlamydia_trachomatis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlamydia_trachomatis");
  }
else if (species == Chlamydophila_abortus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlamydophila_abortus");
  }
else if (species == Chlamydophila_caviae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlamydophila_caviae");
  }
else if (species == Chlamydophila_felis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlamydophila_felis");
  }
else if (species == Chlamydophila_pneumoniae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlamydophila_pneumoniae");
  }
else if (species == Chlamydophila_psittaci)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlamydophila_psittaci");
  }
else if (species == Chlorobaculum_parvum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlorobaculum_parvum");
  }
else if (species == Chlorobaculum_tepidum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlorobaculum_tepidum");
  }
else if (species == Chlorobium_chlorochromatii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlorobium_chlorochromatii");
  }
else if (species == Chlorobium_ferrooxidans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlorobium_ferrooxidans");
  }
else if (species == Chlorobium_limicola)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlorobium_limicola");
  }
else if (species == Chlorobium_phaeobacteroides)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlorobium_phaeobacteroides");
  }
else if (species == Chlorobium_vibrioformis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chlorobium_vibrioformis");
  }
else if (species == Chloroflexus_aggregans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chloroflexus_aggregans");
  }
else if (species == Chloroflexus_aurantiacus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Chloroflexus_aurantiacus");
  }
else if (species == Citrobacter_koseri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Citrobacter_koseri");
  }
else if (species == Citrobacter_rodentium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Citrobacter_rodentium");
  }
else if (species == Citrobacter_youngae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Citrobacter_youngae");
  }
else if (species == Clostridium_acetobutylicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_acetobutylicum");
  }
else if (species == Clostridium_asparagiforme)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_asparagiforme");
  }
else if (species == Clostridium_bartlettii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_bartlettii");
  }
else if (species == Clostridium_beijerinckii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_beijerinckii");
  }
else if (species == Clostridium_bolteae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_bolteae");
  }
else if (species == Clostridium_botulinum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_botulinum");
  }
else if (species == Clostridium_butyricum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_butyricum");
  }
else if (species == Clostridium_carboxidivorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_carboxidivorans");
  }
else if (species == Clostridium_cellulolyticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_cellulolyticum");
  }
else if (species == Clostridium_cellulovorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_cellulovorans");
  }
else if (species == Clostridium_cf)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_cf");
  }
else if (species == Clostridium_difficile)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_difficile");
  }
else if (species == Clostridium_hathewayi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_hathewayi");
  }
else if (species == Clostridium_hiranonis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_hiranonis");
  }
else if (species == Clostridium_hylemonae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_hylemonae");
  }
else if (species == Clostridium_kluyveri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_kluyveri");
  }
else if (species == Clostridium_leptum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_leptum");
  }
else if (species == Clostridium_ljungdahlii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_ljungdahlii");
  }
else if (species == Clostridium_methylpentosum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_methylpentosum");
  }
else if (species == Clostridium_nexile)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_nexile");
  }
else if (species == Clostridium_novyi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_novyi");
  }
else if (species == Clostridium_papyrosolvens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_papyrosolvens");
  }
else if (species == Clostridium_perfringens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_perfringens");
  }
else if (species == Clostridium_phytofermentans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_phytofermentans");
  }
else if (species == Clostridium_saccharolyticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_saccharolyticum");
  }
else if (species == Clostridium_scindens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_scindens");
  }
else if (species == Clostridium_sporogenes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_sporogenes");
  }
else if (species == Clostridium_sticklandii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_sticklandii");
  }
else if (species == Clostridium_symbiosum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_symbiosum");
  }
else if (species == Clostridium_tetani)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_tetani");
  }
else if (species == Clostridium_thermocellum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Clostridium_thermocellum");
  }
else if (species == Collinsella_aerofaciens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Collinsella_aerofaciens");
  }
else if (species == Collinsella_intestinalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Collinsella_intestinalis");
  }
else if (species == Collinsella_stercoris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Collinsella_stercoris");
  }
else if (species == Coprobacillus_bacterium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Coprobacillus_bacterium");
  }
else if (species == Coprococcus_catus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Coprococcus_catus");
  }
else if (species == Coprococcus_comes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Coprococcus_comes");
  }
else if (species == Coprococcus_eutactus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Coprococcus_eutactus");
  }
else if (species == Corynebacterium_accolens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_accolens");
  }
else if (species == Corynebacterium_ammoniagenes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_ammoniagenes");
  }
else if (species == Corynebacterium_amycolatum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_amycolatum");
  }
else if (species == Corynebacterium_aurimucosum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_aurimucosum");
  }
else if (species == Corynebacterium_diphtheriae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_diphtheriae");
  }
else if (species == Corynebacterium_efficiens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_efficiens");
  }
else if (species == Corynebacterium_genitalium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_genitalium");
  }
else if (species == Corynebacterium_glucuronolyticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_glucuronolyticum");
  }
else if (species == Corynebacterium_glutamicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_glutamicum");
  }
else if (species == Corynebacterium_jeikeium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_jeikeium");
  }
else if (species == Corynebacterium_kroppenstedtii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_kroppenstedtii");
  }
else if (species == Corynebacterium_lipophiloflavum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_lipophiloflavum");
  }
else if (species == Corynebacterium_matruchotii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_matruchotii");
  }
else if (species == Corynebacterium_pseudogenitalium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_pseudogenitalium");
  }
else if (species == Corynebacterium_pseudotuberculosis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_pseudotuberculosis");
  }
else if (species == Corynebacterium_resistens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_resistens");
  }
else if (species == Corynebacterium_striatum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_striatum");
  }
else if (species == Corynebacterium_tuberculostearicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_tuberculostearicum");
  }
else if (species == Corynebacterium_urealyticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_urealyticum");
  }
else if (species == Corynebacterium_variabile)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Corynebacterium_variabile");
  }
else if (species == Cronobacter_sakazakii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Cronobacter_sakazakii");
  }
else if (species == Cronobacter_turicensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Cronobacter_turicensis");
  }
else if (species == Cupriavidus_eutropha)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Cupriavidus_eutropha");
  }
else if (species == Cupriavidus_metallidurans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Cupriavidus_metallidurans");
  }
else if (species == Cupriavidus_taiwanensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Cupriavidus_taiwanensis");
  }
else if (species == Dehalococcoides_ethenogenes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Dehalococcoides_ethenogenes");
  }
else if (species == Deinococcus_deserti)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Deinococcus_deserti");
  }
else if (species == Deinococcus_geothermalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Deinococcus_geothermalis");
  }
else if (species == Deinococcus_maricopensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Deinococcus_maricopensis");
  }
else if (species == Deinococcus_proteolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Deinococcus_proteolyticus");
  }
else if (species == Deinococcus_radiodurans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Deinococcus_radiodurans");
  }
else if (species == Desulfotomaculum_acetoxidans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfotomaculum_acetoxidans");
  }
else if (species == Desulfotomaculum_nigrificans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfotomaculum_nigrificans");
  }
else if (species == Desulfotomaculum_reducens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfotomaculum_reducens");
  }
else if (species == Desulfovibrio_aespoeensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfovibrio_aespoeensis");
  }
else if (species == Desulfovibrio_desulfuricans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfovibrio_desulfuricans");
  }
else if (species == Desulfovibrio_fructosovorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfovibrio_fructosovorans");
  }
else if (species == Desulfovibrio_magneticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfovibrio_magneticus");
  }
else if (species == Desulfovibrio_piger)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfovibrio_piger");
  }
else if (species == Desulfovibrio_salexigens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfovibrio_salexigens");
  }
else if (species == Desulfovibrio_vulgaris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfovibrio_vulgaris");
  }
else if (species == Desulfurococcus_kamchatkensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfurococcus_kamchatkensis");
  }
else if (species == Desulfurococcus_mucosus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Desulfurococcus_mucosus");
  }
else if (species == Dialister_invisus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Dialister_invisus");
  }
else if (species == Dialister_microaerophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Dialister_microaerophilus");
  }
else if (species == Dickeya_dadantii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Dickeya_dadantii");
  }
else if (species == Dickeya_zeae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Dickeya_zeae");
  }
else if (species == Dictyoglomus_thermophilum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Dictyoglomus_thermophilum");
  }
else if (species == Dictyoglomus_turgidum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Dictyoglomus_turgidum");
  }
else if (species == Dorea_formicigenerans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Dorea_formicigenerans");
  }
else if (species == Dorea_longicatena)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Dorea_longicatena");
  }
else if (species == Edwardsiella_ictaluri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Edwardsiella_ictaluri");
  }
else if (species == Edwardsiella_tarda)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Edwardsiella_tarda");
  }
else if (species == Eggerthella_lenta)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eggerthella_lenta");
  }
else if (species == Ehrlichia_canis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ehrlichia_canis");
  }
else if (species == Ehrlichia_chaffeensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ehrlichia_chaffeensis");
  }
else if (species == Ehrlichia_ruminantium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ehrlichia_ruminantium");
  }
else if (species == Ensifer_medicae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ensifer_medicae");
  }
else if (species == Ensifer_meliloti)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ensifer_meliloti");
  }
else if (species == Enterobacter_cancerogenus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Enterobacter_cancerogenus");
  }
else if (species == Enterobacter_cloacae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Enterobacter_cloacae");
  }
else if (species == Enterococcus_casseliflavus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Enterococcus_casseliflavus");
  }
else if (species == Enterococcus_faecalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Enterococcus_faecalis");
  }
else if (species == Enterococcus_faecium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Enterococcus_faecium");
  }
else if (species == Enterococcus_gallinarum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Enterococcus_gallinarum");
  }
else if (species == Enterococcus_italicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Enterococcus_italicus");
  }
else if (species == Erwinia_amylovora)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Erwinia_amylovora");
  }
else if (species == Erwinia_billingiae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Erwinia_billingiae");
  }
else if (species == Erwinia_pyrifoliae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Erwinia_pyrifoliae");
  }
else if (species == Erwinia_tasmaniensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Erwinia_tasmaniensis");
  }
else if (species == Erythrobacter_litoralis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Erythrobacter_litoralis");
  }
else if (species == Escherichia_albertii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Escherichia_albertii");
  }
else if (species == Escherichia_coli)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Escherichia_coli");
  }
else if (species == Escherichia_fergusonii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Escherichia_fergusonii");
  }
else if (species == Eubacterium_cellulosolvens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eubacterium_cellulosolvens");
  }
else if (species == Eubacterium_eligens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eubacterium_eligens");
  }
else if (species == Eubacterium_hallii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eubacterium_hallii");
  }
else if (species == Eubacterium_limosum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eubacterium_limosum");
  }
else if (species == Eubacterium_rectale)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eubacterium_rectale");
  }
else if (species == Eubacterium_saburreum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eubacterium_saburreum");
  }
else if (species == Eubacterium_saphenum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eubacterium_saphenum");
  }
else if (species == Eubacterium_siraeum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eubacterium_siraeum");
  }
else if (species == Eubacterium_ventriosum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Eubacterium_ventriosum");
  }
else if (species == Exiguobacterium_sibiricum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Exiguobacterium_sibiricum");
  }
else if (species == Faecalibacterium_cf)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Faecalibacterium_cf");
  }
else if (species == Faecalibacterium_prausnitzii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Faecalibacterium_prausnitzii");
  }
else if (species == Flavobacterium_johnsoniae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Flavobacterium_johnsoniae");
  }
else if (species == Flavobacterium_psychrophilum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Flavobacterium_psychrophilum");
  }
else if (species == Francisella_novicida)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Francisella_novicida");
  }
else if (species == Francisella_philomiragia)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Francisella_philomiragia");
  }
else if (species == Francisella_tularensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Francisella_tularensis");
  }
else if (species == Frankia_alni)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Frankia_alni");
  }
else if (species == Frankia_symbiont)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Frankia_symbiont");
  }
else if (species == Fusobacterium_gonidiaformans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Fusobacterium_gonidiaformans");
  }
else if (species == Fusobacterium_mortiferum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Fusobacterium_mortiferum");
  }
else if (species == Fusobacterium_nucleatum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Fusobacterium_nucleatum");
  }
else if (species == Fusobacterium_periodonticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Fusobacterium_periodonticum");
  }
else if (species == Fusobacterium_ulcerans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Fusobacterium_ulcerans");
  }
else if (species == Fusobacterium_varium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Fusobacterium_varium");
  }
else if (species == Gemella_haemolysans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Gemella_haemolysans");
  }
else if (species == Gemella_moribillum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Gemella_moribillum");
  }
else if (species == Geobacillus_kaustophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Geobacillus_kaustophilus");
  }
else if (species == Geobacillus_thermodenitrificans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Geobacillus_thermodenitrificans");
  }
else if (species == Geobacillus_thermoglucosidasius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Geobacillus_thermoglucosidasius");
  }
else if (species == Geobacter_bemidjiensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Geobacter_bemidjiensis");
  }
else if (species == Geobacter_lovleyi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Geobacter_lovleyi");
  }
else if (species == Geobacter_metallireducens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Geobacter_metallireducens");
  }
else if (species == Geobacter_sulfurreducens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Geobacter_sulfurreducens");
  }
else if (species == Geobacter_uraniumreducens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Geobacter_uraniumreducens");
  }
else if (species == Gluconacetobacter_diazotrophicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Gluconacetobacter_diazotrophicus");
  }
else if (species == Gluconacetobacter_hansenii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Gluconacetobacter_hansenii");
  }
else if (species == Granulicatella_adiacens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Granulicatella_adiacens");
  }
else if (species == Granulicatella_elegans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Granulicatella_elegans");
  }
else if (species == Haemophilus_ducreyi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Haemophilus_ducreyi");
  }
else if (species == Haemophilus_influenzae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Haemophilus_influenzae");
  }
else if (species == Haemophilus_parainfluenzae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Haemophilus_parainfluenzae");
  }
else if (species == Haemophilus_parasuis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Haemophilus_parasuis");
  }
else if (species == Halanaerobium_praevalens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Halanaerobium_praevalens");
  }
else if (species == Halobacterium_salinarum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Halobacterium_salinarum");
  }
else if (species == Helicobacter_acinonychis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_acinonychis");
  }
else if (species == Helicobacter_bilis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_bilis");
  }
else if (species == Helicobacter_canadensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_canadensis");
  }
else if (species == Helicobacter_cinaedi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_cinaedi");
  }
else if (species == Helicobacter_felis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_felis");
  }
else if (species == Helicobacter_hepaticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_hepaticus");
  }
else if (species == Helicobacter_mustelae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_mustelae");
  }
else if (species == Helicobacter_pullorum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_pullorum");
  }
else if (species == Helicobacter_pylori)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_pylori");
  }
else if (species == Helicobacter_suis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_suis");
  }
else if (species == Helicobacter_winghamensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Helicobacter_winghamensis");
  }
else if (species == Idiomarina_baltica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Idiomarina_baltica");
  }
else if (species == Idiomarina_loihiensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Idiomarina_loihiensis");
  }
else if (species == Kingella_denitrificans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Kingella_denitrificans");
  }
else if (species == Kingella_oralis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Kingella_oralis");
  }
else if (species == Klebsiella_pneumoniae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Klebsiella_pneumoniae");
  }
else if (species == Klebsiella_variicola)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Klebsiella_variicola");
  }
else if (species == Labrenzia_aggregata)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Labrenzia_aggregata");
  }
else if (species == Labrenzia_alexandrii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Labrenzia_alexandrii");
  }
else if (species == Lactobacillus_acidophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_acidophilus");
  }
else if (species == Lactobacillus_amylolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_amylolyticus");
  }
else if (species == Lactobacillus_amylovorus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_amylovorus");
  }
else if (species == Lactobacillus_antri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_antri");
  }
else if (species == Lactobacillus_brevis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_brevis");
  }
else if (species == Lactobacillus_buchneri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_buchneri");
  }
else if (species == Lactobacillus_casei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_casei");
  }
else if (species == Lactobacillus_coleohominis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_coleohominis");
  }
else if (species == Lactobacillus_crispatus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_crispatus");
  }
else if (species == Lactobacillus_delbrueckii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_delbrueckii");
  }
else if (species == Lactobacillus_fermentum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_fermentum");
  }
else if (species == Lactobacillus_gasseri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_gasseri");
  }
else if (species == Lactobacillus_helveticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_helveticus");
  }
else if (species == Lactobacillus_hilgardii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_hilgardii");
  }
else if (species == Lactobacillus_iners)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_iners");
  }
else if (species == Lactobacillus_jensenii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_jensenii");
  }
else if (species == Lactobacillus_johnsonii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_johnsonii");
  }
else if (species == Lactobacillus_oris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_oris");
  }
else if (species == Lactobacillus_paracasei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_paracasei");
  }
else if (species == Lactobacillus_plantarum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_plantarum");
  }
else if (species == Lactobacillus_reuteri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_reuteri");
  }
else if (species == Lactobacillus_rhamnosus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_rhamnosus");
  }
else if (species == Lactobacillus_ruminis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_ruminis");
  }
else if (species == Lactobacillus_sakei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_sakei");
  }
else if (species == Lactobacillus_salivarius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_salivarius");
  }
else if (species == Lactobacillus_ultunensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_ultunensis");
  }
else if (species == Lactobacillus_vaginalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lactobacillus_vaginalis");
  }
else if (species == Legionella_drancourtii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Legionella_drancourtii");
  }
else if (species == Legionella_longbeachae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Legionella_longbeachae");
  }
else if (species == Legionella_pneumophila)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Legionella_pneumophila");
  }
else if (species == Leptospira_biflexa)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leptospira_biflexa");
  }
else if (species == Leptospira_borgpetersenii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leptospira_borgpetersenii");
  }
else if (species == Leptospira_interrogans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leptospira_interrogans");
  }
else if (species == Leptotrichia_buccalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leptotrichia_buccalis");
  }
else if (species == Leptotrichia_goodfellowii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leptotrichia_goodfellowii");
  }
else if (species == Leptotrichia_hofstadii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leptotrichia_hofstadii");
  }
else if (species == Leuconostoc_citreum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leuconostoc_citreum");
  }
else if (species == Leuconostoc_gasicomitatum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leuconostoc_gasicomitatum");
  }
else if (species == Leuconostoc_kimchii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leuconostoc_kimchii");
  }
else if (species == Leuconostoc_mesenteroides)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Leuconostoc_mesenteroides");
  }
else if (species == Listeria_grayi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Listeria_grayi");
  }
else if (species == Listeria_innocua)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Listeria_innocua");
  }
else if (species == Listeria_monocytogenes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Listeria_monocytogenes");
  }
else if (species == Listeria_seeligeri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Listeria_seeligeri");
  }
else if (species == Listeria_welshimeri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Listeria_welshimeri");
  }
else if (species == Loktanella_vestfoldensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Loktanella_vestfoldensis");
  }
else if (species == Lysinibacillus_fusiformis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lysinibacillus_fusiformis");
  }
else if (species == Lysinibacillus_sphaericus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Lysinibacillus_sphaericus");
  }
else if (species == Magnetospirillum_magneticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Magnetospirillum_magneticum");
  }
else if (species == Magnetospirillum_magnetotacticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Magnetospirillum_magnetotacticum");
  }
else if (species == Marinobacter_algicola)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Marinobacter_algicola");
  }
else if (species == Marinobacter_aquaeolei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Marinobacter_aquaeolei");
  }
else if (species == Marinobacter_bacterium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Marinobacter_bacterium");
  }
else if (species == Megasphaera_micronuciformis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Megasphaera_micronuciformis");
  }
else if (species == Meiothermus_ruber)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Meiothermus_ruber");
  }
else if (species == Meiothermus_silvanus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Meiothermus_silvanus");
  }
else if (species == Mesorhizobium_ciceri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mesorhizobium_ciceri");
  }
else if (species == Mesorhizobium_loti)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mesorhizobium_loti");
  }
else if (species == Mesorhizobium_opportunistum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mesorhizobium_opportunistum");
  }
else if (species == Methanobrevibacter_ruminantium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanobrevibacter_ruminantium");
  }
else if (species == Methanobrevibacter_smithii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanobrevibacter_smithii");
  }
else if (species == Methanocaldococcus_fervens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanocaldococcus_fervens");
  }
else if (species == Methanocaldococcus_infernus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanocaldococcus_infernus");
  }
else if (species == Methanocaldococcus_jannaschii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanocaldococcus_jannaschii");
  }
else if (species == Methanocaldococcus_vulcanius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanocaldococcus_vulcanius");
  }
else if (species == Methanocella_paludicola)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanocella_paludicola");
  }
else if (species == Methanococcus_aeolicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanococcus_aeolicus");
  }
else if (species == Methanococcus_maripaludis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanococcus_maripaludis");
  }
else if (species == Methanococcus_vannielii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanococcus_vannielii");
  }
else if (species == Methanococcus_voltae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanococcus_voltae");
  }
else if (species == Methanosarcina_acetivorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanosarcina_acetivorans");
  }
else if (species == Methanosarcina_barkeri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanosarcina_barkeri");
  }
else if (species == Methanosarcina_mazei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanosarcina_mazei");
  }
else if (species == Methanothermobacter_marburgensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanothermobacter_marburgensis");
  }
else if (species == Methanothermobacter_thermautotrophicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methanothermobacter_thermautotrophicus");
  }
else if (species == Methylobacterium_chloromethanicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methylobacterium_chloromethanicum");
  }
else if (species == Methylobacterium_extorquens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methylobacterium_extorquens");
  }
else if (species == Methylobacterium_nodulans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methylobacterium_nodulans");
  }
else if (species == Methylobacterium_populi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methylobacterium_populi");
  }
else if (species == Methylobacterium_radiotolerans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methylobacterium_radiotolerans");
  }
else if (species == Methylotenera_mobilis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Methylotenera_mobilis");
  }
else if (species == Micromonospora_aurantiaca)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Micromonospora_aurantiaca");
  }
else if (species == Mobiluncus_curtisii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mobiluncus_curtisii");
  }
else if (species == Mobiluncus_mulieris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mobiluncus_mulieris");
  }
else if (species == Mycobacterium_abscessus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_abscessus");
  }
else if (species == Mycobacterium_avium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_avium");
  }
else if (species == Mycobacterium_bovis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_bovis");
  }
else if (species == Mycobacterium_gilvum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_gilvum");
  }
else if (species == Mycobacterium_intracellulare)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_intracellulare");
  }
else if (species == Mycobacterium_kansasii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_kansasii");
  }
else if (species == Mycobacterium_leprae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_leprae");
  }
else if (species == Mycobacterium_marinum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_marinum");
  }
else if (species == Mycobacterium_parascrofulaceum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_parascrofulaceum");
  }
else if (species == Mycobacterium_smegmatis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_smegmatis");
  }
else if (species == Mycobacterium_tuberculosis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_tuberculosis");
  }
else if (species == Mycobacterium_ulcerans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_ulcerans");
  }
else if (species == Mycobacterium_vanbaalenii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycobacterium_vanbaalenii");
  }
else if (species == Mycoplasma_agalactiae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_agalactiae");
  }
else if (species == Mycoplasma_alligatoris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_alligatoris");
  }
else if (species == Mycoplasma_arthritidis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_arthritidis");
  }
else if (species == Mycoplasma_bovis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_bovis");
  }
else if (species == Mycoplasma_capricolum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_capricolum");
  }
else if (species == Mycoplasma_conjunctivae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_conjunctivae");
  }
else if (species == Mycoplasma_crocodyli)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_crocodyli");
  }
else if (species == Mycoplasma_fermentans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_fermentans");
  }
else if (species == Mycoplasma_gallisepticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_gallisepticum");
  }
else if (species == Mycoplasma_genitalium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_genitalium");
  }
else if (species == Mycoplasma_haemofelis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_haemofelis");
  }
else if (species == Mycoplasma_hominis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_hominis");
  }
else if (species == Mycoplasma_hyopneumoniae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_hyopneumoniae");
  }
else if (species == Mycoplasma_hyorhinis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_hyorhinis");
  }
else if (species == Mycoplasma_leachii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_leachii");
  }
else if (species == Mycoplasma_mobile)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_mobile");
  }
else if (species == Mycoplasma_mycoides)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_mycoides");
  }
else if (species == Mycoplasma_penetrans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_penetrans");
  }
else if (species == Mycoplasma_pneumoniae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_pneumoniae");
  }
else if (species == Mycoplasma_pulmonis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_pulmonis");
  }
else if (species == Mycoplasma_synoviae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Mycoplasma_synoviae");
  }
else if (species == Neisseria_cinerea)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_cinerea");
  }
else if (species == Neisseria_elongata)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_elongata");
  }
else if (species == Neisseria_flavescens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_flavescens");
  }
else if (species == Neisseria_gonorrhoeae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_gonorrhoeae");
  }
else if (species == Neisseria_lactamica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_lactamica");
  }
else if (species == Neisseria_meningitidis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_meningitidis");
  }
else if (species == Neisseria_mucosa)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_mucosa");
  }
else if (species == Neisseria_polysaccharea)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_polysaccharea");
  }
else if (species == Neisseria_sicca)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_sicca");
  }
else if (species == Neisseria_subflava)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neisseria_subflava");
  }
else if (species == Neorickettsia_risticii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neorickettsia_risticii");
  }
else if (species == Neorickettsia_sennetsu)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Neorickettsia_sennetsu");
  }
else if (species == Nitrobacter_hamburgensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Nitrobacter_hamburgensis");
  }
else if (species == Nitrobacter_winogradskyi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Nitrobacter_winogradskyi");
  }
else if (species == Nitrosococcus_halophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Nitrosococcus_halophilus");
  }
else if (species == Nitrosococcus_oceani)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Nitrosococcus_oceani");
  }
else if (species == Nitrosococcus_watsoni)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Nitrosococcus_watsoni");
  }
else if (species == Nitrosomonas_europaea)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Nitrosomonas_europaea");
  }
else if (species == Nitrosomonas_eutropha)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Nitrosomonas_eutropha");
  }
else if (species == Nostoc_punctiforme)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Nostoc_punctiforme");
  }
else if (species == Oceanicola_batsensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Oceanicola_batsensis");
  }
else if (species == Oceanicola_granulosus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Oceanicola_granulosus");
  }
else if (species == Ochrobactrum_anthropi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ochrobactrum_anthropi");
  }
else if (species == Ochrobactrum_intermedium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ochrobactrum_intermedium");
  }
else if (species == Oribacterium_sinus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Oribacterium_sinus");
  }
else if (species == Paenibacillus_curdlanolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Paenibacillus_curdlanolyticus");
  }
else if (species == Paenibacillus_larvae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Paenibacillus_larvae");
  }
else if (species == Paenibacillus_polymyxa)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Paenibacillus_polymyxa");
  }
else if (species == Paenibacillus_vortex)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Paenibacillus_vortex");
  }
else if (species == Pantoea_ananatis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pantoea_ananatis");
  }
else if (species == Pantoea_vagans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pantoea_vagans");
  }
else if (species == Parabacteroides_distasonis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Parabacteroides_distasonis");
  }
else if (species == Parabacteroides_johnsonii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Parabacteroides_johnsonii");
  }
else if (species == Parabacteroides_merdae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Parabacteroides_merdae");
  }
else if (species == Pasteurella_dagmatis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pasteurella_dagmatis");
  }
else if (species == Pasteurella_multocida)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pasteurella_multocida");
  }
else if (species == Pectobacterium_carotovora)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pectobacterium_carotovora");
  }
else if (species == Pectobacterium_carotovorum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pectobacterium_carotovorum");
  }
else if (species == Pectobacterium_wasabiae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pectobacterium_wasabiae");
  }
else if (species == Pediococcus_acidilactici)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pediococcus_acidilactici");
  }
else if (species == Pediococcus_pentosaceus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pediococcus_pentosaceus");
  }
else if (species == Pedobacter_heparinus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pedobacter_heparinus");
  }
else if (species == Pedobacter_saltans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pedobacter_saltans");
  }
else if (species == Pelobacter_carbinolicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pelobacter_carbinolicus");
  }
else if (species == Pelobacter_propionicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pelobacter_propionicus");
  }
else if (species == Pelodictyon_clathratiforme)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pelodictyon_clathratiforme");
  }
else if (species == Pelodictyon_luteolum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pelodictyon_luteolum");
  }
else if (species == Peptoniphilus_duerdenii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Peptoniphilus_duerdenii");
  }
else if (species == Peptoniphilus_harei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Peptoniphilus_harei");
  }
else if (species == Peptoniphilus_lacrimalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Peptoniphilus_lacrimalis");
  }
else if (species == Peptostreptococcus_anaerobius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Peptostreptococcus_anaerobius");
  }
else if (species == Peptostreptococcus_stomatis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Peptostreptococcus_stomatis");
  }
else if (species == Photobacterium_angustum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Photobacterium_angustum");
  }
else if (species == Photobacterium_damselae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Photobacterium_damselae");
  }
else if (species == Photobacterium_profundum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Photobacterium_profundum");
  }
else if (species == Photorhabdus_asymbiotica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Photorhabdus_asymbiotica");
  }
else if (species == Photorhabdus_luminescens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Photorhabdus_luminescens");
  }
else if (species == Planctomyces_brasiliensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Planctomyces_brasiliensis");
  }
else if (species == Planctomyces_limnophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Planctomyces_limnophilus");
  }
else if (species == Planctomyces_maris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Planctomyces_maris");
  }
else if (species == Polaribacter_irgensii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Polaribacter_irgensii");
  }
else if (species == Polaromonas_naphthalenivorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Polaromonas_naphthalenivorans");
  }
else if (species == Polynucleobacter_necessarius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Polynucleobacter_necessarius");
  }
else if (species == Porphyromonas_asaccharolytica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Porphyromonas_asaccharolytica");
  }
else if (species == Porphyromonas_endodontalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Porphyromonas_endodontalis");
  }
else if (species == Porphyromonas_gingivalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Porphyromonas_gingivalis");
  }
else if (species == Porphyromonas_uenonis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Porphyromonas_uenonis");
  }
else if (species == Prevotella_amnii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_amnii");
  }
else if (species == Prevotella_bergensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_bergensis");
  }
else if (species == Prevotella_bivia)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_bivia");
  }
else if (species == Prevotella_bryantii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_bryantii");
  }
else if (species == Prevotella_buccae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_buccae");
  }
else if (species == Prevotella_buccalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_buccalis");
  }
else if (species == Prevotella_copri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_copri");
  }
else if (species == Prevotella_disiens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_disiens");
  }
else if (species == Prevotella_marshii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_marshii");
  }
else if (species == Prevotella_melaninogenica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_melaninogenica");
  }
else if (species == Prevotella_multiformis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_multiformis");
  }
else if (species == Prevotella_oralis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_oralis");
  }
else if (species == Prevotella_oris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_oris");
  }
else if (species == Prevotella_ruminicola)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_ruminicola");
  }
else if (species == Prevotella_salivae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_salivae");
  }
else if (species == Prevotella_tannerae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_tannerae");
  }
else if (species == Prevotella_timonensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_timonensis");
  }
else if (species == Prevotella_veroralis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Prevotella_veroralis");
  }
else if (species == Propionibacterium_acnes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Propionibacterium_acnes");
  }
else if (species == Propionibacterium_freudenreichii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Propionibacterium_freudenreichii");
  }
else if (species == Proteus_mirabilis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Proteus_mirabilis");
  }
else if (species == Proteus_penneri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Proteus_penneri");
  }
else if (species == Providencia_alcalifaciens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Providencia_alcalifaciens");
  }
else if (species == Providencia_rettgeri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Providencia_rettgeri");
  }
else if (species == Providencia_rustigianii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Providencia_rustigianii");
  }
else if (species == Providencia_stuartii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Providencia_stuartii");
  }
else if (species == Pseudoalteromonas_atlantica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudoalteromonas_atlantica");
  }
else if (species == Pseudoalteromonas_haloplanktis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudoalteromonas_haloplanktis");
  }
else if (species == Pseudoalteromonas_tunicata)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudoalteromonas_tunicata");
  }
else if (species == Pseudomonas_aeruginosa)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudomonas_aeruginosa");
  }
else if (species == Pseudomonas_entomophila)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudomonas_entomophila");
  }
else if (species == Pseudomonas_fluorescens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudomonas_fluorescens");
  }
else if (species == Pseudomonas_mendocina)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudomonas_mendocina");
  }
else if (species == Pseudomonas_putida)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudomonas_putida");
  }
else if (species == Pseudomonas_savastanoi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudomonas_savastanoi");
  }
else if (species == Pseudomonas_stutzeri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudomonas_stutzeri");
  }
else if (species == Pseudomonas_syringae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pseudomonas_syringae");
  }
else if (species == Psychrobacter_arcticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Psychrobacter_arcticus");
  }
else if (species == Psychrobacter_cryohalolentis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Psychrobacter_cryohalolentis");
  }
else if (species == Psychromonas_ingrahamii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Psychromonas_ingrahamii");
  }
else if (species == Pyrobaculum_aerophilum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pyrobaculum_aerophilum");
  }
else if (species == Pyrobaculum_arsenaticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pyrobaculum_arsenaticum");
  }
else if (species == Pyrobaculum_calidifontis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pyrobaculum_calidifontis");
  }
else if (species == Pyrobaculum_islandicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pyrobaculum_islandicum");
  }
else if (species == Pyrococcus_abyssi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pyrococcus_abyssi");
  }
else if (species == Pyrococcus_furiosus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pyrococcus_furiosus");
  }
else if (species == Pyrococcus_horikoshii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Pyrococcus_horikoshii");
  }
else if (species == Ralstonia_pickettii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ralstonia_pickettii");
  }
else if (species == Ralstonia_solanacearum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ralstonia_solanacearum");
  }
else if (species == Rhizobium_etli)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rhizobium_etli");
  }
else if (species == Rhizobium_leguminosarum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rhizobium_leguminosarum");
  }
else if (species == Rhizobium_radiobacter)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rhizobium_radiobacter");
  }
else if (species == Rhodobacter_capsulatus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rhodobacter_capsulatus");
  }
else if (species == Rhodobacter_sphaeroides)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rhodobacter_sphaeroides");
  }
else if (species == Rhodococcus_equi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rhodococcus_equi");
  }
else if (species == Rhodococcus_erythropolis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rhodococcus_erythropolis");
  }
else if (species == Rhodococcus_opacus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rhodococcus_opacus");
  }
else if (species == Rickettsia_africae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_africae");
  }
else if (species == Rickettsia_akari)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_akari");
  }
else if (species == Rickettsia_bellii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_bellii");
  }
else if (species == Rickettsia_canadensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_canadensis");
  }
else if (species == Rickettsia_conorii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_conorii");
  }
else if (species == Rickettsia_endosymbiont)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_endosymbiont");
  }
else if (species == Rickettsia_felis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_felis");
  }
else if (species == Rickettsia_massiliae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_massiliae");
  }
else if (species == Rickettsia_peacockii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_peacockii");
  }
else if (species == Rickettsia_prowazekii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_prowazekii");
  }
else if (species == Rickettsia_rickettsii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_rickettsii");
  }
else if (species == Rickettsia_sibirica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_sibirica");
  }
else if (species == Rickettsia_typhi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rickettsia_typhi");
  }
else if (species == Roseburia_intestinalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Roseburia_intestinalis");
  }
else if (species == Roseburia_inulinivorans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Roseburia_inulinivorans");
  }
else if (species == Roseiflexus_castenholzii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Roseiflexus_castenholzii");
  }
else if (species == Roseobacter_denitrificans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Roseobacter_denitrificans");
  }
else if (species == Roseobacter_litoralis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Roseobacter_litoralis");
  }
else if (species == Roseovarius_nubinhibens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Roseovarius_nubinhibens");
  }
else if (species == Rothia_dentocariosa)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rothia_dentocariosa");
  }
else if (species == Rothia_mucilaginosa)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Rothia_mucilaginosa");
  }
else if (species == Ruegeria_lacuscaerulensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ruegeria_lacuscaerulensis");
  }
else if (species == Ruegeria_pomeroyi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ruegeria_pomeroyi");
  }
else if (species == Ruminococcus_albus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ruminococcus_albus");
  }
else if (species == Ruminococcus_bromii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ruminococcus_bromii");
  }
else if (species == Ruminococcus_flavefaciens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ruminococcus_flavefaciens");
  }
else if (species == Ruminococcus_gnavus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ruminococcus_gnavus");
  }
else if (species == Ruminococcus_lactaris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ruminococcus_lactaris");
  }
else if (species == Ruminococcus_obeum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ruminococcus_obeum");
  }
else if (species == Ruminococcus_torques)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ruminococcus_torques");
  }
else if (species == Salinispora_arenicola)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Salinispora_arenicola");
  }
else if (species == Salinispora_tropica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Salinispora_tropica");
  }
else if (species == Salmonella_enterica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Salmonella_enterica");
  }
else if (species == Salmonella_typhimurium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Salmonella_typhimurium");
  }
else if (species == Segniliparus_rotundus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Segniliparus_rotundus");
  }
else if (species == Segniliparus_rugosus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Segniliparus_rugosus");
  }
else if (species == Selenomonas_artemidis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Selenomonas_artemidis");
  }
else if (species == Selenomonas_flueggei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Selenomonas_flueggei");
  }
else if (species == Selenomonas_noxia)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Selenomonas_noxia");
  }
else if (species == Selenomonas_sputigena)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Selenomonas_sputigena");
  }
else if (species == Serratia_odorifera)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Serratia_odorifera");
  }
else if (species == Serratia_proteamaculans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Serratia_proteamaculans");
  }
else if (species == Shewanella_amazonensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_amazonensis");
  }
else if (species == Shewanella_baltica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_baltica");
  }
else if (species == Shewanella_benthica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_benthica");
  }
else if (species == Shewanella_denitrificans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_denitrificans");
  }
else if (species == Shewanella_frigidimarina)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_frigidimarina");
  }
else if (species == Shewanella_halifaxensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_halifaxensis");
  }
else if (species == Shewanella_loihica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_loihica");
  }
else if (species == Shewanella_oneidensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_oneidensis");
  }
else if (species == Shewanella_pealeana)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_pealeana");
  }
else if (species == Shewanella_piezotolerans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_piezotolerans");
  }
else if (species == Shewanella_putrefaciens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_putrefaciens");
  }
else if (species == Shewanella_sediminis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_sediminis");
  }
else if (species == Shewanella_violacea)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_violacea");
  }
else if (species == Shewanella_woodyi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shewanella_woodyi");
  }
else if (species == Shigella_boydii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shigella_boydii");
  }
else if (species == Shigella_dysenteriae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shigella_dysenteriae");
  }
else if (species == Shigella_flexneri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shigella_flexneri");
  }
else if (species == Shigella_sonnei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shigella_sonnei");
  }
else if (species == Slackia_exigua)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Slackia_exigua");
  }
else if (species == Slackia_heliotrinireducens)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Slackia_heliotrinireducens");
  }
else if (species == Sphingobium_chlorophenolicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sphingobium_chlorophenolicum");
  }
else if (species == Sphingobium_japonicum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sphingobium_japonicum");
  }
else if (species == Sphingomonas_wittichii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sphingomonas_wittichii");
  }
else if (species == Spirochaeta_smaragdinae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Spirochaeta_smaragdinae");
  }
else if (species == Spirochaeta_thermophila)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Spirochaeta_thermophila");
  }
else if (species == Saureus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Saureus");
  }
else if (species == Scapitis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Scapitis");
  }
else if (species == Scarnosus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Scarnosus");
  }
else if (species == Sepidermidis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sepidermidis");
  }
else if (species == Shaemolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shaemolyticus");
  }
else if (species == Shominis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Shominis");
  }
else if (species == Slugdunensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Slugdunensis");
  }
else if (species == Spseudintermedius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Spseudintermedius");
  }
else if (species == Ssaprophyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ssaprophyticus");
  }
else if (species == Swarneri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Swarneri");
  }
else if (species == Staphylothermus_hellenicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Staphylothermus_hellenicus");
  }
else if (species == Staphylothermus_marinus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Staphylothermus_marinus");
  }
else if (species == Stenotrophomonas_maltophilia)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Stenotrophomonas_maltophilia");
  }
else if (species == Streptococcus_agalactiae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_agalactiae");
  }
else if (species == Streptococcus_anginosus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_anginosus");
  }
else if (species == Streptococcus_australis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_australis");
  }
else if (species == Streptococcus_bovis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_bovis");
  }
else if (species == Streptococcus_cristatus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_cristatus");
  }
else if (species == Streptococcus_downei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_downei");
  }
else if (species == Streptococcus_dysgalactiae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_dysgalactiae");
  }
else if (species == Streptococcus_equi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_equi");
  }
else if (species == Streptococcus_equinus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_equinus");
  }
else if (species == Streptococcus_gallolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_gallolyticus");
  }
else if (species == Streptococcus_gordonii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_gordonii");
  }
else if (species == Streptococcus_infantarius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_infantarius");
  }
else if (species == Streptococcus_infantis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_infantis");
  }
else if (species == Streptococcus_mitis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_mitis");
  }
else if (species == Streptococcus_mutans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_mutans");
  }
else if (species == Streptococcus_oralis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_oralis");
  }
else if (species == Streptococcus_parasanguinis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_parasanguinis");
  }
else if (species == Streptococcus_peroris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_peroris");
  }
else if (species == Streptococcus_pneumoniae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_pneumoniae");
  }
else if (species == Streptococcus_pseudoporcinus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_pseudoporcinus");
  }
else if (species == Streptococcus_pyogenes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_pyogenes");
  }
else if (species == Streptococcus_salivarius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_salivarius");
  }
else if (species == Streptococcus_sanguinis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_sanguinis");
  }
else if (species == Streptococcus_suis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_suis");
  }
else if (species == Streptococcus_thermophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_thermophilus");
  }
else if (species == Streptococcus_uberis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_uberis");
  }
else if (species == Streptococcus_vestibularis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptococcus_vestibularis");
  }
else if (species == Streptomyces_albus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_albus");
  }
else if (species == Streptomyces_avermitilis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_avermitilis");
  }
else if (species == Streptomyces_bingchenggensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_bingchenggensis");
  }
else if (species == Streptomyces_clavuligerus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_clavuligerus");
  }
else if (species == Streptomyces_coelicolor)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_coelicolor");
  }
else if (species == Streptomyces_ghanaensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_ghanaensis");
  }
else if (species == Streptomyces_griseoflavus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_griseoflavus");
  }
else if (species == Streptomyces_griseus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_griseus");
  }
else if (species == Streptomyces_hygroscopicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_hygroscopicus");
  }
else if (species == Streptomyces_lividans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_lividans");
  }
else if (species == Streptomyces_pristinaespiralis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_pristinaespiralis");
  }
else if (species == Streptomyces_roseosporus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_roseosporus");
  }
else if (species == Streptomyces_scabiei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_scabiei");
  }
else if (species == Streptomyces_sviceus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_sviceus");
  }
else if (species == Streptomyces_violaceusniger)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_violaceusniger");
  }
else if (species == Streptomyces_viridochromogenes)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Streptomyces_viridochromogenes");
  }
else if (species == Sulfolobus_acidocaldarius)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sulfolobus_acidocaldarius");
  }
else if (species == Sulfolobus_islandicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sulfolobus_islandicus");
  }
else if (species == Sulfolobus_solfataricus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sulfolobus_solfataricus");
  }
else if (species == Sulfolobus_tokodaii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sulfolobus_tokodaii");
  }
else if (species == Sulfurihydrogenibium_azorense)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sulfurihydrogenibium_azorense");
  }
else if (species == Sulfurihydrogenibium_yellowstonense)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sulfurihydrogenibium_yellowstonense");
  }
else if (species == Sulfurimonas_autotrophica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sulfurimonas_autotrophica");
  }
else if (species == Sulfurimonas_denitrificans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Sulfurimonas_denitrificans");
  }
else if (species == Synechococcus_elongatus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Synechococcus_elongatus");
  }
else if (species == Thermaerobacter_marianensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermaerobacter_marianensis");
  }
else if (species == Thermaerobacter_subterraneus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermaerobacter_subterraneus");
  }
else if (species == Thermoanaerobacter_brockii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoanaerobacter_brockii");
  }
else if (species == Thermoanaerobacter_ethanolicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoanaerobacter_ethanolicus");
  }
else if (species == Thermoanaerobacter_italicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoanaerobacter_italicus");
  }
else if (species == Thermoanaerobacterium_thermosaccharolyticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoanaerobacterium_thermosaccharolyticum");
  }
else if (species == Thermoanaerobacterium_xylanolyticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoanaerobacterium_xylanolyticum");
  }
else if (species == Thermoanaerobacter_mathranii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoanaerobacter_mathranii");
  }
else if (species == Thermoanaerobacter_pseudethanolicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoanaerobacter_pseudethanolicus");
  }
else if (species == Thermoanaerobacter_wiegelii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoanaerobacter_wiegelii");
  }
else if (species == Thermococcus_barophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermococcus_barophilus");
  }
else if (species == Thermococcus_gammatolerans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermococcus_gammatolerans");
  }
else if (species == Thermococcus_kodakaraensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermococcus_kodakaraensis");
  }
else if (species == Thermococcus_onnurineus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermococcus_onnurineus");
  }
else if (species == Thermococcus_sibiricus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermococcus_sibiricus");
  }
else if (species == Thermoplasma_acidophilum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoplasma_acidophilum");
  }
else if (species == Thermoplasma_volcanium)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermoplasma_volcanium");
  }
else if (species == Thermosipho_africanus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermosipho_africanus");
  }
else if (species == Thermosipho_melanesiensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermosipho_melanesiensis");
  }
else if (species == Thermotoga_lettingae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermotoga_lettingae");
  }
else if (species == Thermotoga_maritima)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermotoga_maritima");
  }
else if (species == Thermotoga_naphthophila)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermotoga_naphthophila");
  }
else if (species == Thermotoga_neapolitana)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermotoga_neapolitana");
  }
else if (species == Thermotoga_petrophila)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermotoga_petrophila");
  }
else if (species == Thermus_aquaticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermus_aquaticus");
  }
else if (species == Thermus_scotoductus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermus_scotoductus");
  }
else if (species == Thermus_thermophilus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Thermus_thermophilus");
  }
else if (species == Treponema_denticola)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Treponema_denticola");
  }
else if (species == Treponema_pallidum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Treponema_pallidum");
  }
else if (species == Treponema_phagedenis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Treponema_phagedenis");
  }
else if (species == Treponema_vincentii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Treponema_vincentii");
  }
else if (species == Ureaplasma_parvum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ureaplasma_parvum");
  }
else if (species == Ureaplasma_urealyticum)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Ureaplasma_urealyticum");
  }
else if (species == Veillonella_atypica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Veillonella_atypica");
  }
else if (species == Veillonella_dispar)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Veillonella_dispar");
  }
else if (species == Veillonella_parvula)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Veillonella_parvula");
  }
else if (species == Vibrio_alginolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_alginolyticus");
  }
else if (species == Vibrio_brasiliensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_brasiliensis");
  }
else if (species == Vibrio_campbellii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_campbellii");
  }
else if (species == Vibrio_caribbenthicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_caribbenthicus");
  }
else if (species == Vibrio_cholerae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_cholerae");
  }
else if (species == Vibrio_coralliilyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_coralliilyticus");
  }
else if (species == Vibrio_furnissii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_furnissii");
  }
else if (species == Vibrio_harveyi)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_harveyi");
  }
else if (species == Vibrio_metschnikovii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_metschnikovii");
  }
else if (species == Vibrio_mimicus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_mimicus");
  }
else if (species == Vibrio_orientalis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_orientalis");
  }
else if (species == Vibrio_parahaemolyticus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_parahaemolyticus");
  }
else if (species == Vibrio_shilonii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_shilonii");
  }
else if (species == Vibrio_sinaloensis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_sinaloensis");
  }
else if (species == Vibrio_splendidus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_splendidus");
  }
else if (species == Vibrio_vulnificus)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vibrio_vulnificus");
  }
else if (species == Vulcanisaeta_distributa)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vulcanisaeta_distributa");
  }
else if (species == Vulcanisaeta_moutnovskia)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Vulcanisaeta_moutnovskia");
  }
else if (species == Wolbachia_endosymbiont)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Wolbachia_endosymbiont");
  }
else if (species == Wolbachia_pipientis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Wolbachia_pipientis");
  }
else if (species == Xanthomonas_albilineans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Xanthomonas_albilineans");
  }
else if (species == Xanthomonas_axonopodis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Xanthomonas_axonopodis");
  }
else if (species == Xanthomonas_campestris)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Xanthomonas_campestris");
  }
else if (species == Xanthomonas_fuscans)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Xanthomonas_fuscans");
  }
else if (species == Xanthomonas_oryzae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Xanthomonas_oryzae");
  }
else if (species == Xenorhabdus_bovienii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Xenorhabdus_bovienii");
  }
else if (species == Xenorhabdus_nematophila)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Xenorhabdus_nematophila");
  }
else if (species == Yersinia_aldovae)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_aldovae");
  }
else if (species == Yersinia_bercovieri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_bercovieri");
  }
else if (species == Yersinia_enterocolitica)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_enterocolitica");
  }
else if (species == Yersinia_frederiksenii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_frederiksenii");
  }
else if (species == Yersinia_intermedia)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_intermedia");
  }
else if (species == Yersinia_kristensenii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_kristensenii");
  }
else if (species == Yersinia_mollaretii)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_mollaretii");
  }
else if (species == Yersinia_pestis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_pestis");
  }
else if (species == Yersinia_pseudotuberculosis)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_pseudotuberculosis");
  }
else if (species == Yersinia_rohdei)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_rohdei");
  }
else if (species == Yersinia_ruckeri)
  {
    strbuf_reset(sbuf);
    strbuf_append_str(sbuf, "Yersinia_ruckeri");
  }      
  else
    {
      die("Coding error - I would expect the compiler to prevent assigning a bad enum value - we get %d\n", species);
    }
  
}



void load_all_phylo_group_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
  panel_file_paths[CoagPos] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[CoagPos], "data/staph/species/Saureus.fasta" );
  panel_file_paths[CoagNeg] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[CoagNeg], "data/staph/species/coag_neg.fasta" );  
}

// void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
// {
//   panel_file_paths[Saureus] = strbuf_create(install_dir->buff);
//   strbuf_append_str(panel_file_paths[Saureus], "data/staph/species/Saureus.fasta" );
//   panel_file_paths[Sepidermidis] = strbuf_create(install_dir->buff);
//   strbuf_append_str(panel_file_paths[Sepidermidis], "data/staph/species/Sepidermidis.fasta" );
//   panel_file_paths[Shaemolyticus] = strbuf_create(install_dir->buff);
//   strbuf_append_str(panel_file_paths[Shaemolyticus], "data/staph/species/Shaemolyticus.fasta" );
//   panel_file_paths[Sother] = strbuf_create(install_dir->buff);
//   strbuf_append_str(panel_file_paths[Sother], "data/staph/species/Sother.fasta" );
// }

void load_all_species_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
  panel_file_paths[Sother] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sother], "data/staph/species/Sother.fasta" );  
  panel_file_paths[Achromobacter_piechaudii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Achromobacter_piechaudii], "data/meta_species/species/Achromobacter_piechaudii.fa");
  panel_file_paths[Achromobacter_xylosoxidans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Achromobacter_xylosoxidans], "data/meta_species/species/Achromobacter_xylosoxidans.fa");
  panel_file_paths[Acidaminococcus_fermentans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acidaminococcus_fermentans], "data/meta_species/species/Acidaminococcus_fermentans.fa");
  panel_file_paths[Acidithiobacillus_caldus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acidithiobacillus_caldus], "data/meta_species/species/Acidithiobacillus_caldus.fa");
  panel_file_paths[Acidithiobacillus_ferrooxidans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acidithiobacillus_ferrooxidans], "data/meta_species/species/Acidithiobacillus_ferrooxidans.fa");
  panel_file_paths[Acidovorax_avenae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acidovorax_avenae], "data/meta_species/species/Acidovorax_avenae.fa");
  panel_file_paths[Acidovorax_delafieldii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acidovorax_delafieldii], "data/meta_species/species/Acidovorax_delafieldii.fa");
  panel_file_paths[Acinetobacter_baumannii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acinetobacter_baumannii], "data/meta_species/species/Acinetobacter_baumannii.fa");
  panel_file_paths[Acinetobacter_calcoaceticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acinetobacter_calcoaceticus], "data/meta_species/species/Acinetobacter_calcoaceticus.fa");
  panel_file_paths[Acinetobacter_haemolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acinetobacter_haemolyticus], "data/meta_species/species/Acinetobacter_haemolyticus.fa");
  panel_file_paths[Acinetobacter_johnsonii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acinetobacter_johnsonii], "data/meta_species/species/Acinetobacter_johnsonii.fa");
  panel_file_paths[Acinetobacter_junii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acinetobacter_junii], "data/meta_species/species/Acinetobacter_junii.fa");
  panel_file_paths[Acinetobacter_lwoffii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acinetobacter_lwoffii], "data/meta_species/species/Acinetobacter_lwoffii.fa");
  panel_file_paths[Acinetobacter_radioresistens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Acinetobacter_radioresistens], "data/meta_species/species/Acinetobacter_radioresistens.fa");
  panel_file_paths[Actinobacillus_minor] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Actinobacillus_minor], "data/meta_species/species/Actinobacillus_minor.fa");
  panel_file_paths[Actinobacillus_pleuropneumoniae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Actinobacillus_pleuropneumoniae], "data/meta_species/species/Actinobacillus_pleuropneumoniae.fa");
  panel_file_paths[Actinobacillus_succinogenes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Actinobacillus_succinogenes], "data/meta_species/species/Actinobacillus_succinogenes.fa");
  panel_file_paths[Actinobacillus_ureae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Actinobacillus_ureae], "data/meta_species/species/Actinobacillus_ureae.fa");
  panel_file_paths[Actinomyces_coleocanis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Actinomyces_coleocanis], "data/meta_species/species/Actinomyces_coleocanis.fa");
  panel_file_paths[Actinomyces_odontolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Actinomyces_odontolyticus], "data/meta_species/species/Actinomyces_odontolyticus.fa");
  panel_file_paths[Actinomyces_oris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Actinomyces_oris], "data/meta_species/species/Actinomyces_oris.fa");
  panel_file_paths[Actinomyces_urogenitalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Actinomyces_urogenitalis], "data/meta_species/species/Actinomyces_urogenitalis.fa");
  panel_file_paths[Actinomyces_viscosus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Actinomyces_viscosus], "data/meta_species/species/Actinomyces_viscosus.fa");
  panel_file_paths[Aeromonas_hydrophila] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Aeromonas_hydrophila], "data/meta_species/species/Aeromonas_hydrophila.fa");
  panel_file_paths[Aeromonas_salmonicida] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Aeromonas_salmonicida], "data/meta_species/species/Aeromonas_salmonicida.fa");
  panel_file_paths[Aggregatibacter_actinomycetemcomitans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Aggregatibacter_actinomycetemcomitans], "data/meta_species/species/Aggregatibacter_actinomycetemcomitans.fa");
  panel_file_paths[Aggregatibacter_aphrophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Aggregatibacter_aphrophilus], "data/meta_species/species/Aggregatibacter_aphrophilus.fa");
  panel_file_paths[Aggregatibacter_segnis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Aggregatibacter_segnis], "data/meta_species/species/Aggregatibacter_segnis.fa");
  panel_file_paths[Agrobacterium_tumefaciens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Agrobacterium_tumefaciens], "data/meta_species/species/Agrobacterium_tumefaciens.fa");
  panel_file_paths[Agrobacterium_vitis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Agrobacterium_vitis], "data/meta_species/species/Agrobacterium_vitis.fa");
  panel_file_paths[Alcanivorax_borkumensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Alcanivorax_borkumensis], "data/meta_species/species/Alcanivorax_borkumensis.fa");
  panel_file_paths[Alicyclobacillus_acidocaldarius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Alicyclobacillus_acidocaldarius], "data/meta_species/species/Alicyclobacillus_acidocaldarius.fa");
  panel_file_paths[Alicyclobacillus_Alicyclobacillus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Alicyclobacillus_Alicyclobacillus], "data/meta_species/species/Alicyclobacillus_Alicyclobacillus.fa");
  panel_file_paths[Aliivibrio_fischeri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Aliivibrio_fischeri], "data/meta_species/species/Aliivibrio_fischeri.fa");
  panel_file_paths[Aliivibrio_salmonicida] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Aliivibrio_salmonicida], "data/meta_species/species/Aliivibrio_salmonicida.fa");
  panel_file_paths[Alistipes_putredinis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Alistipes_putredinis], "data/meta_species/species/Alistipes_putredinis.fa");
  panel_file_paths[Alistipes_shahii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Alistipes_shahii], "data/meta_species/species/Alistipes_shahii.fa");
  panel_file_paths[Alkaliphilus_metalliredigens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Alkaliphilus_metalliredigens], "data/meta_species/species/Alkaliphilus_metalliredigens.fa");
  panel_file_paths[Alkaliphilus_oremlandii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Alkaliphilus_oremlandii], "data/meta_species/species/Alkaliphilus_oremlandii.fa");
  panel_file_paths[Anabaena_azollae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anabaena_azollae], "data/meta_species/species/Anabaena_azollae.fa");
  panel_file_paths[Anabaena_variabilis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anabaena_variabilis], "data/meta_species/species/Anabaena_variabilis.fa");
  panel_file_paths[Anaerococcus_hydrogenalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaerococcus_hydrogenalis], "data/meta_species/species/Anaerococcus_hydrogenalis.fa");
  panel_file_paths[Anaerococcus_lactolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaerococcus_lactolyticus], "data/meta_species/species/Anaerococcus_lactolyticus.fa");
  panel_file_paths[Anaerococcus_prevotii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaerococcus_prevotii], "data/meta_species/species/Anaerococcus_prevotii.fa");
  panel_file_paths[Anaerococcus_tetradius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaerococcus_tetradius], "data/meta_species/species/Anaerococcus_tetradius.fa");
  panel_file_paths[Anaerococcus_vaginalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaerococcus_vaginalis], "data/meta_species/species/Anaerococcus_vaginalis.fa");
  panel_file_paths[Anaeromyxobacter_dehalogenans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaeromyxobacter_dehalogenans], "data/meta_species/species/Anaeromyxobacter_dehalogenans.fa");
  panel_file_paths[Anaerostipes_caccae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaerostipes_caccae], "data/meta_species/species/Anaerostipes_caccae.fa");
  panel_file_paths[Anaplasma_centrale] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaplasma_centrale], "data/meta_species/species/Anaplasma_centrale.fa");
  panel_file_paths[Anaplasma_marginale] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaplasma_marginale], "data/meta_species/species/Anaplasma_marginale.fa");
  panel_file_paths[Anaplasma_phagocytophilum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Anaplasma_phagocytophilum], "data/meta_species/species/Anaplasma_phagocytophilum.fa");
  panel_file_paths[Archaeoglobus_fulgidus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Archaeoglobus_fulgidus], "data/meta_species/species/Archaeoglobus_fulgidus.fa");
  panel_file_paths[Archaeoglobus_profundus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Archaeoglobus_profundus], "data/meta_species/species/Archaeoglobus_profundus.fa");
  panel_file_paths[Arcobacter_butzleri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Arcobacter_butzleri], "data/meta_species/species/Arcobacter_butzleri.fa");
  panel_file_paths[Arcobacter_nitrofigilis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Arcobacter_nitrofigilis], "data/meta_species/species/Arcobacter_nitrofigilis.fa");
  panel_file_paths[Arthrobacter_arilaitensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Arthrobacter_arilaitensis], "data/meta_species/species/Arthrobacter_arilaitensis.fa");
  panel_file_paths[Arthrobacter_aurescens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Arthrobacter_aurescens], "data/meta_species/species/Arthrobacter_aurescens.fa");
  panel_file_paths[Arthrobacter_chlorophenolicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Arthrobacter_chlorophenolicus], "data/meta_species/species/Arthrobacter_chlorophenolicus.fa");
  panel_file_paths[Arthrobacter_phenanthrenivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Arthrobacter_phenanthrenivorans], "data/meta_species/species/Arthrobacter_phenanthrenivorans.fa");
  panel_file_paths[Arthrospira_maxima] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Arthrospira_maxima], "data/meta_species/species/Arthrospira_maxima.fa");
  panel_file_paths[Arthrospira_platensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Arthrospira_platensis], "data/meta_species/species/Arthrospira_platensis.fa");
  panel_file_paths[Atopobium_parvulum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Atopobium_parvulum], "data/meta_species/species/Atopobium_parvulum.fa");
  panel_file_paths[Atopobium_rimae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Atopobium_rimae], "data/meta_species/species/Atopobium_rimae.fa");
  panel_file_paths[Atopobium_vaginae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Atopobium_vaginae], "data/meta_species/species/Atopobium_vaginae.fa");
  panel_file_paths[Bacillus_amyloliquefaciens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_amyloliquefaciens], "data/meta_species/species/Bacillus_amyloliquefaciens.fa");
  panel_file_paths[Bacillus_anthracis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_anthracis], "data/meta_species/species/Bacillus_anthracis.fa");
  panel_file_paths[Bacillus_atrophaeus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_atrophaeus], "data/meta_species/species/Bacillus_atrophaeus.fa");
  panel_file_paths[Bacillus_cellulosilyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_cellulosilyticus], "data/meta_species/species/Bacillus_cellulosilyticus.fa");
  panel_file_paths[Bacillus_cereus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_cereus], "data/meta_species/species/Bacillus_cereus.fa");
  panel_file_paths[Bacillus_clausii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_clausii], "data/meta_species/species/Bacillus_clausii.fa");
  panel_file_paths[Bacillus_coagulans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_coagulans], "data/meta_species/species/Bacillus_coagulans.fa");
  panel_file_paths[Bacillus_coahuilensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_coahuilensis], "data/meta_species/species/Bacillus_coahuilensis.fa");
  panel_file_paths[Bacillus_halodurans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_halodurans], "data/meta_species/species/Bacillus_halodurans.fa");
  panel_file_paths[Bacillus_licheniformis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_licheniformis], "data/meta_species/species/Bacillus_licheniformis.fa");
  panel_file_paths[Bacillus_megaterium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_megaterium], "data/meta_species/species/Bacillus_megaterium.fa");
  panel_file_paths[Bacillus_mycoides] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_mycoides], "data/meta_species/species/Bacillus_mycoides.fa");
  panel_file_paths[Bacillus_pseudofirmus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_pseudofirmus], "data/meta_species/species/Bacillus_pseudofirmus.fa");
  panel_file_paths[Bacillus_pseudomycoides] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_pseudomycoides], "data/meta_species/species/Bacillus_pseudomycoides.fa");
  panel_file_paths[Bacillus_pumilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_pumilus], "data/meta_species/species/Bacillus_pumilus.fa");
  panel_file_paths[Bacillus_selenitireducens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_selenitireducens], "data/meta_species/species/Bacillus_selenitireducens.fa");
  panel_file_paths[Bacillus_subtilis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_subtilis], "data/meta_species/species/Bacillus_subtilis.fa");
  panel_file_paths[Bacillus_thuringiensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_thuringiensis], "data/meta_species/species/Bacillus_thuringiensis.fa");
  panel_file_paths[Bacillus_weihenstephanensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacillus_weihenstephanensis], "data/meta_species/species/Bacillus_weihenstephanensis.fa");
  panel_file_paths[Bacteroides_caccae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_caccae], "data/meta_species/species/Bacteroides_caccae.fa");
  panel_file_paths[Bacteroides_cellulosilyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_cellulosilyticus], "data/meta_species/species/Bacteroides_cellulosilyticus.fa");
  panel_file_paths[Bacteroides_coprocola] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_coprocola], "data/meta_species/species/Bacteroides_coprocola.fa");
  panel_file_paths[Bacteroides_coprophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_coprophilus], "data/meta_species/species/Bacteroides_coprophilus.fa");
  panel_file_paths[Bacteroides_dorei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_dorei], "data/meta_species/species/Bacteroides_dorei.fa");
  panel_file_paths[Bacteroides_eggerthii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_eggerthii], "data/meta_species/species/Bacteroides_eggerthii.fa");
  panel_file_paths[Bacteroides_finegoldii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_finegoldii], "data/meta_species/species/Bacteroides_finegoldii.fa");
  panel_file_paths[Bacteroides_fragilis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_fragilis], "data/meta_species/species/Bacteroides_fragilis.fa");
  panel_file_paths[Bacteroides_helcogenes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_helcogenes], "data/meta_species/species/Bacteroides_helcogenes.fa");
  panel_file_paths[Bacteroides_intestinalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_intestinalis], "data/meta_species/species/Bacteroides_intestinalis.fa");
  panel_file_paths[Bacteroides_ovatus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_ovatus], "data/meta_species/species/Bacteroides_ovatus.fa");
  panel_file_paths[Bacteroides_pectinophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_pectinophilus], "data/meta_species/species/Bacteroides_pectinophilus.fa");
  panel_file_paths[Bacteroides_plebeius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_plebeius], "data/meta_species/species/Bacteroides_plebeius.fa");
  panel_file_paths[Bacteroides_salanitronis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_salanitronis], "data/meta_species/species/Bacteroides_salanitronis.fa");
  panel_file_paths[Bacteroides_stercoris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_stercoris], "data/meta_species/species/Bacteroides_stercoris.fa");
  panel_file_paths[Bacteroides_thetaiotaomicron] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_thetaiotaomicron], "data/meta_species/species/Bacteroides_thetaiotaomicron.fa");
  panel_file_paths[Bacteroides_uniformis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_uniformis], "data/meta_species/species/Bacteroides_uniformis.fa");
  panel_file_paths[Bacteroides_vulgatus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_vulgatus], "data/meta_species/species/Bacteroides_vulgatus.fa");
  panel_file_paths[Bacteroides_xylanisolvens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bacteroides_xylanisolvens], "data/meta_species/species/Bacteroides_xylanisolvens.fa");
  panel_file_paths[Bartonella_bacilliformis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bartonella_bacilliformis], "data/meta_species/species/Bartonella_bacilliformis.fa");
  panel_file_paths[Bartonella_clarridgeiae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bartonella_clarridgeiae], "data/meta_species/species/Bartonella_clarridgeiae.fa");
  panel_file_paths[Bartonella_grahamii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bartonella_grahamii], "data/meta_species/species/Bartonella_grahamii.fa");
  panel_file_paths[Bartonella_henselae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bartonella_henselae], "data/meta_species/species/Bartonella_henselae.fa");
  panel_file_paths[Bartonella_quintana] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bartonella_quintana], "data/meta_species/species/Bartonella_quintana.fa");
  panel_file_paths[Bartonella_tribocorum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bartonella_tribocorum], "data/meta_species/species/Bartonella_tribocorum.fa");
  panel_file_paths[Bifidobacterium_adolescentis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_adolescentis], "data/meta_species/species/Bifidobacterium_adolescentis.fa");
  panel_file_paths[Bifidobacterium_angulatum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_angulatum], "data/meta_species/species/Bifidobacterium_angulatum.fa");
  panel_file_paths[Bifidobacterium_animalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_animalis], "data/meta_species/species/Bifidobacterium_animalis.fa");
  panel_file_paths[Bifidobacterium_bifidum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_bifidum], "data/meta_species/species/Bifidobacterium_bifidum.fa");
  panel_file_paths[Bifidobacterium_breve] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_breve], "data/meta_species/species/Bifidobacterium_breve.fa");
  panel_file_paths[Bifidobacterium_catenulatum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_catenulatum], "data/meta_species/species/Bifidobacterium_catenulatum.fa");
  panel_file_paths[Bifidobacterium_dentium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_dentium], "data/meta_species/species/Bifidobacterium_dentium.fa");
  panel_file_paths[Bifidobacterium_gallicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_gallicum], "data/meta_species/species/Bifidobacterium_gallicum.fa");
  panel_file_paths[Bifidobacterium_longum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_longum], "data/meta_species/species/Bifidobacterium_longum.fa");
  panel_file_paths[Bifidobacterium_pseudocatenulatum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bifidobacterium_pseudocatenulatum], "data/meta_species/species/Bifidobacterium_pseudocatenulatum.fa");
  panel_file_paths[Blautia_hansenii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Blautia_hansenii], "data/meta_species/species/Blautia_hansenii.fa");
  panel_file_paths[Blautia_hydrogenotrophica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Blautia_hydrogenotrophica], "data/meta_species/species/Blautia_hydrogenotrophica.fa");
  panel_file_paths[Bordetella_avium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bordetella_avium], "data/meta_species/species/Bordetella_avium.fa");
  panel_file_paths[Bordetella_bronchiseptica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bordetella_bronchiseptica], "data/meta_species/species/Bordetella_bronchiseptica.fa");
  panel_file_paths[Bordetella_parapertussis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bordetella_parapertussis], "data/meta_species/species/Bordetella_parapertussis.fa");
  panel_file_paths[Bordetella_pertussis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bordetella_pertussis], "data/meta_species/species/Bordetella_pertussis.fa");
  panel_file_paths[Bordetella_petrii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bordetella_petrii], "data/meta_species/species/Bordetella_petrii.fa");
  panel_file_paths[Borrelia_afzelii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Borrelia_afzelii], "data/meta_species/species/Borrelia_afzelii.fa");
  panel_file_paths[Borrelia_burgdorferi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Borrelia_burgdorferi], "data/meta_species/species/Borrelia_burgdorferi.fa");
  panel_file_paths[Borrelia_duttonii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Borrelia_duttonii], "data/meta_species/species/Borrelia_duttonii.fa");
  panel_file_paths[Borrelia_garinii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Borrelia_garinii], "data/meta_species/species/Borrelia_garinii.fa");
  panel_file_paths[Borrelia_hermsii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Borrelia_hermsii], "data/meta_species/species/Borrelia_hermsii.fa");
  panel_file_paths[Borrelia_recurrentis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Borrelia_recurrentis], "data/meta_species/species/Borrelia_recurrentis.fa");
  panel_file_paths[Borrelia_spielmanii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Borrelia_spielmanii], "data/meta_species/species/Borrelia_spielmanii.fa");
  panel_file_paths[Borrelia_turicatae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Borrelia_turicatae], "data/meta_species/species/Borrelia_turicatae.fa");
  panel_file_paths[Borrelia_valaisiana] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Borrelia_valaisiana], "data/meta_species/species/Borrelia_valaisiana.fa");
  panel_file_paths[Brachyspira_hyodysenteriae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brachyspira_hyodysenteriae], "data/meta_species/species/Brachyspira_hyodysenteriae.fa");
  panel_file_paths[Brachyspira_murdochii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brachyspira_murdochii], "data/meta_species/species/Brachyspira_murdochii.fa");
  panel_file_paths[Brachyspira_pilosicoli] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brachyspira_pilosicoli], "data/meta_species/species/Brachyspira_pilosicoli.fa");
  panel_file_paths[Bradyrhizobium_japonicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Bradyrhizobium_japonicum], "data/meta_species/species/Bradyrhizobium_japonicum.fa");
  panel_file_paths[Brevibacterium_linens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brevibacterium_linens], "data/meta_species/species/Brevibacterium_linens.fa");
  panel_file_paths[Brevibacterium_mcbrellneri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brevibacterium_mcbrellneri], "data/meta_species/species/Brevibacterium_mcbrellneri.fa");
  panel_file_paths[Brevundimonas_subvibrioides] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brevundimonas_subvibrioides], "data/meta_species/species/Brevundimonas_subvibrioides.fa");
  panel_file_paths[Brucella_abortus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brucella_abortus], "data/meta_species/species/Brucella_abortus.fa");
  panel_file_paths[Brucella_canis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brucella_canis], "data/meta_species/species/Brucella_canis.fa");
  panel_file_paths[Brucella_ceti] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brucella_ceti], "data/meta_species/species/Brucella_ceti.fa");
  panel_file_paths[Brucella_melitensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brucella_melitensis], "data/meta_species/species/Brucella_melitensis.fa");
  panel_file_paths[Brucella_microti] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brucella_microti], "data/meta_species/species/Brucella_microti.fa");
  panel_file_paths[Brucella_neotomae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brucella_neotomae], "data/meta_species/species/Brucella_neotomae.fa");
  panel_file_paths[Brucella_ovis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brucella_ovis], "data/meta_species/species/Brucella_ovis.fa");
  panel_file_paths[Brucella_pinnipedialis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brucella_pinnipedialis], "data/meta_species/species/Brucella_pinnipedialis.fa");
  panel_file_paths[Brucella_suis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Brucella_suis], "data/meta_species/species/Brucella_suis.fa");
  panel_file_paths[Burkholderia_ambifaria] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_ambifaria], "data/meta_species/species/Burkholderia_ambifaria.fa");
  panel_file_paths[Burkholderia_cenocepacia] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_cenocepacia], "data/meta_species/species/Burkholderia_cenocepacia.fa");
  panel_file_paths[Burkholderia_cepacia] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_cepacia], "data/meta_species/species/Burkholderia_cepacia.fa");
  panel_file_paths[Burkholderia_dolosa] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_dolosa], "data/meta_species/species/Burkholderia_dolosa.fa");
  panel_file_paths[Burkholderia_glumae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_glumae], "data/meta_species/species/Burkholderia_glumae.fa");
  panel_file_paths[Burkholderia_graminis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_graminis], "data/meta_species/species/Burkholderia_graminis.fa");
  panel_file_paths[Burkholderia_mallei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_mallei], "data/meta_species/species/Burkholderia_mallei.fa");
  panel_file_paths[Burkholderia_multivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_multivorans], "data/meta_species/species/Burkholderia_multivorans.fa");
  panel_file_paths[Burkholderia_oklahomensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_oklahomensis], "data/meta_species/species/Burkholderia_oklahomensis.fa");
  panel_file_paths[Burkholderia_phymatum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_phymatum], "data/meta_species/species/Burkholderia_phymatum.fa");
  panel_file_paths[Burkholderia_phytofirmans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_phytofirmans], "data/meta_species/species/Burkholderia_phytofirmans.fa");
  panel_file_paths[Burkholderia_pseudomallei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_pseudomallei], "data/meta_species/species/Burkholderia_pseudomallei.fa");
  panel_file_paths[Burkholderia_rhizoxinica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_rhizoxinica], "data/meta_species/species/Burkholderia_rhizoxinica.fa");
  panel_file_paths[Burkholderia_thailandensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_thailandensis], "data/meta_species/species/Burkholderia_thailandensis.fa");
  panel_file_paths[Burkholderia_ubonensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_ubonensis], "data/meta_species/species/Burkholderia_ubonensis.fa");
  panel_file_paths[Burkholderia_vietnamiensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_vietnamiensis], "data/meta_species/species/Burkholderia_vietnamiensis.fa");
  panel_file_paths[Burkholderia_xenovorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Burkholderia_xenovorans], "data/meta_species/species/Burkholderia_xenovorans.fa");
  panel_file_paths[Butyrivibrio_crossotus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Butyrivibrio_crossotus], "data/meta_species/species/Butyrivibrio_crossotus.fa");
  panel_file_paths[Butyrivibrio_fibrisolvens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Butyrivibrio_fibrisolvens], "data/meta_species/species/Butyrivibrio_fibrisolvens.fa");
  panel_file_paths[Butyrivibrio_proteoclasticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Butyrivibrio_proteoclasticus], "data/meta_species/species/Butyrivibrio_proteoclasticus.fa");
  panel_file_paths[Caldanaerobacter_pacificum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldanaerobacter_pacificum], "data/meta_species/species/Caldanaerobacter_pacificum.fa");
  panel_file_paths[Caldanaerobacter_tengcongensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldanaerobacter_tengcongensis], "data/meta_species/species/Caldanaerobacter_tengcongensis.fa");
  panel_file_paths[Caldicellulosiruptor_bescii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldicellulosiruptor_bescii], "data/meta_species/species/Caldicellulosiruptor_bescii.fa");
  panel_file_paths[Caldicellulosiruptor_hydrothermalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldicellulosiruptor_hydrothermalis], "data/meta_species/species/Caldicellulosiruptor_hydrothermalis.fa");
  panel_file_paths[Caldicellulosiruptor_kristjanssonii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldicellulosiruptor_kristjanssonii], "data/meta_species/species/Caldicellulosiruptor_kristjanssonii.fa");
  panel_file_paths[Caldicellulosiruptor_kronotskyensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldicellulosiruptor_kronotskyensis], "data/meta_species/species/Caldicellulosiruptor_kronotskyensis.fa");
  panel_file_paths[Caldicellulosiruptor_lactoaceticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldicellulosiruptor_lactoaceticus], "data/meta_species/species/Caldicellulosiruptor_lactoaceticus.fa");
  panel_file_paths[Caldicellulosiruptor_obsidiansis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldicellulosiruptor_obsidiansis], "data/meta_species/species/Caldicellulosiruptor_obsidiansis.fa");
  panel_file_paths[Caldicellulosiruptor_owensensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldicellulosiruptor_owensensis], "data/meta_species/species/Caldicellulosiruptor_owensensis.fa");
  panel_file_paths[Caldicellulosiruptor_saccharolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caldicellulosiruptor_saccharolyticus], "data/meta_species/species/Caldicellulosiruptor_saccharolyticus.fa");
  panel_file_paths[Campylobacter_coli] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_coli], "data/meta_species/species/Campylobacter_coli.fa");
  panel_file_paths[Campylobacter_concisus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_concisus], "data/meta_species/species/Campylobacter_concisus.fa");
  panel_file_paths[Campylobacter_curvus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_curvus], "data/meta_species/species/Campylobacter_curvus.fa");
  panel_file_paths[Campylobacter_fetus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_fetus], "data/meta_species/species/Campylobacter_fetus.fa");
  panel_file_paths[Campylobacter_gracilis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_gracilis], "data/meta_species/species/Campylobacter_gracilis.fa");
  panel_file_paths[Campylobacter_hominis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_hominis], "data/meta_species/species/Campylobacter_hominis.fa");
  panel_file_paths[Campylobacter_jejuni] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_jejuni], "data/meta_species/species/Campylobacter_jejuni.fa");
  panel_file_paths[Campylobacter_lari] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_lari], "data/meta_species/species/Campylobacter_lari.fa");
  panel_file_paths[Campylobacter_rectus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_rectus], "data/meta_species/species/Campylobacter_rectus.fa");
  panel_file_paths[Campylobacter_showae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_showae], "data/meta_species/species/Campylobacter_showae.fa");
  panel_file_paths[Campylobacter_upsaliensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Campylobacter_upsaliensis], "data/meta_species/species/Campylobacter_upsaliensis.fa");
  panel_file_paths[Candidatus_Blochmannia_floridanus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Candidatus_Blochmannia_floridanus], "data/meta_species/species/Candidatus_Blochmannia_floridanus.fa");
  panel_file_paths[Candidatus_Blochmannia_pennsylvanicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Candidatus_Blochmannia_pennsylvanicus], "data/meta_species/species/Candidatus_Blochmannia_pennsylvanicus.fa");
  panel_file_paths[Candidatus_Pelagibacter_ubique] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Candidatus_Pelagibacter_ubique], "data/meta_species/species/Candidatus_Pelagibacter_ubique.fa");
  panel_file_paths[Candidatus_Phytoplasma_yellows] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Candidatus_Phytoplasma_yellows], "data/meta_species/species/Candidatus_Phytoplasma_yellows.fa");
  panel_file_paths[Candidatus_Sulcia_muelleri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Candidatus_Sulcia_muelleri], "data/meta_species/species/Candidatus_Sulcia_muelleri.fa");
  panel_file_paths[Capnocytophaga_gingivalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Capnocytophaga_gingivalis], "data/meta_species/species/Capnocytophaga_gingivalis.fa");
  panel_file_paths[Capnocytophaga_ochracea] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Capnocytophaga_ochracea], "data/meta_species/species/Capnocytophaga_ochracea.fa");
  panel_file_paths[Capnocytophaga_sputigena] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Capnocytophaga_sputigena], "data/meta_species/species/Capnocytophaga_sputigena.fa");
  panel_file_paths[Caulobacter_crescentus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caulobacter_crescentus], "data/meta_species/species/Caulobacter_crescentus.fa");
  panel_file_paths[Caulobacter_segnis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Caulobacter_segnis], "data/meta_species/species/Caulobacter_segnis.fa");
  panel_file_paths[Cellulophaga_algicola] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Cellulophaga_algicola], "data/meta_species/species/Cellulophaga_algicola.fa");
  panel_file_paths[Cellulophaga_lytica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Cellulophaga_lytica], "data/meta_species/species/Cellulophaga_lytica.fa");
  panel_file_paths[Chlamydia_muridarum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlamydia_muridarum], "data/meta_species/species/Chlamydia_muridarum.fa");
  panel_file_paths[Chlamydia_trachomatis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlamydia_trachomatis], "data/meta_species/species/Chlamydia_trachomatis.fa");
  panel_file_paths[Chlamydophila_abortus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlamydophila_abortus], "data/meta_species/species/Chlamydophila_abortus.fa");
  panel_file_paths[Chlamydophila_caviae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlamydophila_caviae], "data/meta_species/species/Chlamydophila_caviae.fa");
  panel_file_paths[Chlamydophila_felis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlamydophila_felis], "data/meta_species/species/Chlamydophila_felis.fa");
  panel_file_paths[Chlamydophila_pneumoniae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlamydophila_pneumoniae], "data/meta_species/species/Chlamydophila_pneumoniae.fa");
  panel_file_paths[Chlamydophila_psittaci] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlamydophila_psittaci], "data/meta_species/species/Chlamydophila_psittaci.fa");
  panel_file_paths[Chlorobaculum_parvum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlorobaculum_parvum], "data/meta_species/species/Chlorobaculum_parvum.fa");
  panel_file_paths[Chlorobaculum_tepidum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlorobaculum_tepidum], "data/meta_species/species/Chlorobaculum_tepidum.fa");
  panel_file_paths[Chlorobium_chlorochromatii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlorobium_chlorochromatii], "data/meta_species/species/Chlorobium_chlorochromatii.fa");
  panel_file_paths[Chlorobium_ferrooxidans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlorobium_ferrooxidans], "data/meta_species/species/Chlorobium_ferrooxidans.fa");
  panel_file_paths[Chlorobium_limicola] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlorobium_limicola], "data/meta_species/species/Chlorobium_limicola.fa");
  panel_file_paths[Chlorobium_phaeobacteroides] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlorobium_phaeobacteroides], "data/meta_species/species/Chlorobium_phaeobacteroides.fa");
  panel_file_paths[Chlorobium_vibrioformis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chlorobium_vibrioformis], "data/meta_species/species/Chlorobium_vibrioformis.fa");
  panel_file_paths[Chloroflexus_aggregans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chloroflexus_aggregans], "data/meta_species/species/Chloroflexus_aggregans.fa");
  panel_file_paths[Chloroflexus_aurantiacus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Chloroflexus_aurantiacus], "data/meta_species/species/Chloroflexus_aurantiacus.fa");
  panel_file_paths[Citrobacter_koseri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Citrobacter_koseri], "data/meta_species/species/Citrobacter_koseri.fa");
  panel_file_paths[Citrobacter_rodentium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Citrobacter_rodentium], "data/meta_species/species/Citrobacter_rodentium.fa");
  panel_file_paths[Citrobacter_youngae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Citrobacter_youngae], "data/meta_species/species/Citrobacter_youngae.fa");
  panel_file_paths[Clostridium_acetobutylicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_acetobutylicum], "data/meta_species/species/Clostridium_acetobutylicum.fa");
  panel_file_paths[Clostridium_asparagiforme] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_asparagiforme], "data/meta_species/species/Clostridium_asparagiforme.fa");
  panel_file_paths[Clostridium_bartlettii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_bartlettii], "data/meta_species/species/Clostridium_bartlettii.fa");
  panel_file_paths[Clostridium_beijerinckii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_beijerinckii], "data/meta_species/species/Clostridium_beijerinckii.fa");
  panel_file_paths[Clostridium_bolteae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_bolteae], "data/meta_species/species/Clostridium_bolteae.fa");
  panel_file_paths[Clostridium_botulinum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_botulinum], "data/meta_species/species/Clostridium_botulinum.fa");
  panel_file_paths[Clostridium_butyricum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_butyricum], "data/meta_species/species/Clostridium_butyricum.fa");
  panel_file_paths[Clostridium_carboxidivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_carboxidivorans], "data/meta_species/species/Clostridium_carboxidivorans.fa");
  panel_file_paths[Clostridium_cellulolyticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_cellulolyticum], "data/meta_species/species/Clostridium_cellulolyticum.fa");
  panel_file_paths[Clostridium_cellulovorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_cellulovorans], "data/meta_species/species/Clostridium_cellulovorans.fa");
  panel_file_paths[Clostridium_cf] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_cf], "data/meta_species/species/Clostridium_cf.fa");
  panel_file_paths[Clostridium_difficile] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_difficile], "data/meta_species/species/Clostridium_difficile.fa");
  panel_file_paths[Clostridium_hathewayi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_hathewayi], "data/meta_species/species/Clostridium_hathewayi.fa");
  panel_file_paths[Clostridium_hiranonis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_hiranonis], "data/meta_species/species/Clostridium_hiranonis.fa");
  panel_file_paths[Clostridium_hylemonae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_hylemonae], "data/meta_species/species/Clostridium_hylemonae.fa");
  panel_file_paths[Clostridium_kluyveri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_kluyveri], "data/meta_species/species/Clostridium_kluyveri.fa");
  panel_file_paths[Clostridium_leptum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_leptum], "data/meta_species/species/Clostridium_leptum.fa");
  panel_file_paths[Clostridium_ljungdahlii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_ljungdahlii], "data/meta_species/species/Clostridium_ljungdahlii.fa");
  panel_file_paths[Clostridium_methylpentosum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_methylpentosum], "data/meta_species/species/Clostridium_methylpentosum.fa");
  panel_file_paths[Clostridium_nexile] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_nexile], "data/meta_species/species/Clostridium_nexile.fa");
  panel_file_paths[Clostridium_novyi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_novyi], "data/meta_species/species/Clostridium_novyi.fa");
  panel_file_paths[Clostridium_papyrosolvens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_papyrosolvens], "data/meta_species/species/Clostridium_papyrosolvens.fa");
  panel_file_paths[Clostridium_perfringens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_perfringens], "data/meta_species/species/Clostridium_perfringens.fa");
  panel_file_paths[Clostridium_phytofermentans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_phytofermentans], "data/meta_species/species/Clostridium_phytofermentans.fa");
  panel_file_paths[Clostridium_saccharolyticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_saccharolyticum], "data/meta_species/species/Clostridium_saccharolyticum.fa");
  panel_file_paths[Clostridium_scindens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_scindens], "data/meta_species/species/Clostridium_scindens.fa");
  panel_file_paths[Clostridium_sporogenes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_sporogenes], "data/meta_species/species/Clostridium_sporogenes.fa");
  panel_file_paths[Clostridium_sticklandii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_sticklandii], "data/meta_species/species/Clostridium_sticklandii.fa");
  panel_file_paths[Clostridium_symbiosum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_symbiosum], "data/meta_species/species/Clostridium_symbiosum.fa");
  panel_file_paths[Clostridium_tetani] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_tetani], "data/meta_species/species/Clostridium_tetani.fa");
  panel_file_paths[Clostridium_thermocellum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Clostridium_thermocellum], "data/meta_species/species/Clostridium_thermocellum.fa");
  panel_file_paths[Collinsella_aerofaciens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Collinsella_aerofaciens], "data/meta_species/species/Collinsella_aerofaciens.fa");
  panel_file_paths[Collinsella_intestinalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Collinsella_intestinalis], "data/meta_species/species/Collinsella_intestinalis.fa");
  panel_file_paths[Collinsella_stercoris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Collinsella_stercoris], "data/meta_species/species/Collinsella_stercoris.fa");
  panel_file_paths[Coprobacillus_bacterium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Coprobacillus_bacterium], "data/meta_species/species/Coprobacillus_bacterium.fa");
  panel_file_paths[Coprococcus_catus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Coprococcus_catus], "data/meta_species/species/Coprococcus_catus.fa");
  panel_file_paths[Coprococcus_comes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Coprococcus_comes], "data/meta_species/species/Coprococcus_comes.fa");
  panel_file_paths[Coprococcus_eutactus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Coprococcus_eutactus], "data/meta_species/species/Coprococcus_eutactus.fa");
  panel_file_paths[Corynebacterium_accolens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_accolens], "data/meta_species/species/Corynebacterium_accolens.fa");
  panel_file_paths[Corynebacterium_ammoniagenes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_ammoniagenes], "data/meta_species/species/Corynebacterium_ammoniagenes.fa");
  panel_file_paths[Corynebacterium_amycolatum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_amycolatum], "data/meta_species/species/Corynebacterium_amycolatum.fa");
  panel_file_paths[Corynebacterium_aurimucosum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_aurimucosum], "data/meta_species/species/Corynebacterium_aurimucosum.fa");
  panel_file_paths[Corynebacterium_diphtheriae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_diphtheriae], "data/meta_species/species/Corynebacterium_diphtheriae.fa");
  panel_file_paths[Corynebacterium_efficiens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_efficiens], "data/meta_species/species/Corynebacterium_efficiens.fa");
  panel_file_paths[Corynebacterium_genitalium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_genitalium], "data/meta_species/species/Corynebacterium_genitalium.fa");
  panel_file_paths[Corynebacterium_glucuronolyticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_glucuronolyticum], "data/meta_species/species/Corynebacterium_glucuronolyticum.fa");
  panel_file_paths[Corynebacterium_glutamicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_glutamicum], "data/meta_species/species/Corynebacterium_glutamicum.fa");
  panel_file_paths[Corynebacterium_jeikeium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_jeikeium], "data/meta_species/species/Corynebacterium_jeikeium.fa");
  panel_file_paths[Corynebacterium_kroppenstedtii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_kroppenstedtii], "data/meta_species/species/Corynebacterium_kroppenstedtii.fa");
  panel_file_paths[Corynebacterium_lipophiloflavum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_lipophiloflavum], "data/meta_species/species/Corynebacterium_lipophiloflavum.fa");
  panel_file_paths[Corynebacterium_matruchotii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_matruchotii], "data/meta_species/species/Corynebacterium_matruchotii.fa");
  panel_file_paths[Corynebacterium_pseudogenitalium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_pseudogenitalium], "data/meta_species/species/Corynebacterium_pseudogenitalium.fa");
  panel_file_paths[Corynebacterium_pseudotuberculosis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_pseudotuberculosis], "data/meta_species/species/Corynebacterium_pseudotuberculosis.fa");
  panel_file_paths[Corynebacterium_resistens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_resistens], "data/meta_species/species/Corynebacterium_resistens.fa");
  panel_file_paths[Corynebacterium_striatum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_striatum], "data/meta_species/species/Corynebacterium_striatum.fa");
  panel_file_paths[Corynebacterium_tuberculostearicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_tuberculostearicum], "data/meta_species/species/Corynebacterium_tuberculostearicum.fa");
  panel_file_paths[Corynebacterium_urealyticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_urealyticum], "data/meta_species/species/Corynebacterium_urealyticum.fa");
  panel_file_paths[Corynebacterium_variabile] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Corynebacterium_variabile], "data/meta_species/species/Corynebacterium_variabile.fa");
  panel_file_paths[Cronobacter_sakazakii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Cronobacter_sakazakii], "data/meta_species/species/Cronobacter_sakazakii.fa");
  panel_file_paths[Cronobacter_turicensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Cronobacter_turicensis], "data/meta_species/species/Cronobacter_turicensis.fa");
  panel_file_paths[Cupriavidus_eutropha] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Cupriavidus_eutropha], "data/meta_species/species/Cupriavidus_eutropha.fa");
  panel_file_paths[Cupriavidus_metallidurans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Cupriavidus_metallidurans], "data/meta_species/species/Cupriavidus_metallidurans.fa");
  panel_file_paths[Cupriavidus_taiwanensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Cupriavidus_taiwanensis], "data/meta_species/species/Cupriavidus_taiwanensis.fa");
  panel_file_paths[Dehalococcoides_ethenogenes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Dehalococcoides_ethenogenes], "data/meta_species/species/Dehalococcoides_ethenogenes.fa");
  panel_file_paths[Deinococcus_deserti] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Deinococcus_deserti], "data/meta_species/species/Deinococcus_deserti.fa");
  panel_file_paths[Deinococcus_geothermalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Deinococcus_geothermalis], "data/meta_species/species/Deinococcus_geothermalis.fa");
  panel_file_paths[Deinococcus_maricopensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Deinococcus_maricopensis], "data/meta_species/species/Deinococcus_maricopensis.fa");
  panel_file_paths[Deinococcus_proteolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Deinococcus_proteolyticus], "data/meta_species/species/Deinococcus_proteolyticus.fa");
  panel_file_paths[Deinococcus_radiodurans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Deinococcus_radiodurans], "data/meta_species/species/Deinococcus_radiodurans.fa");
  panel_file_paths[Desulfotomaculum_acetoxidans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfotomaculum_acetoxidans], "data/meta_species/species/Desulfotomaculum_acetoxidans.fa");
  panel_file_paths[Desulfotomaculum_nigrificans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfotomaculum_nigrificans], "data/meta_species/species/Desulfotomaculum_nigrificans.fa");
  panel_file_paths[Desulfotomaculum_reducens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfotomaculum_reducens], "data/meta_species/species/Desulfotomaculum_reducens.fa");
  panel_file_paths[Desulfovibrio_aespoeensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfovibrio_aespoeensis], "data/meta_species/species/Desulfovibrio_aespoeensis.fa");
  panel_file_paths[Desulfovibrio_desulfuricans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfovibrio_desulfuricans], "data/meta_species/species/Desulfovibrio_desulfuricans.fa");
  panel_file_paths[Desulfovibrio_fructosovorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfovibrio_fructosovorans], "data/meta_species/species/Desulfovibrio_fructosovorans.fa");
  panel_file_paths[Desulfovibrio_magneticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfovibrio_magneticus], "data/meta_species/species/Desulfovibrio_magneticus.fa");
  panel_file_paths[Desulfovibrio_piger] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfovibrio_piger], "data/meta_species/species/Desulfovibrio_piger.fa");
  panel_file_paths[Desulfovibrio_salexigens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfovibrio_salexigens], "data/meta_species/species/Desulfovibrio_salexigens.fa");
  panel_file_paths[Desulfovibrio_vulgaris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfovibrio_vulgaris], "data/meta_species/species/Desulfovibrio_vulgaris.fa");
  panel_file_paths[Desulfurococcus_kamchatkensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfurococcus_kamchatkensis], "data/meta_species/species/Desulfurococcus_kamchatkensis.fa");
  panel_file_paths[Desulfurococcus_mucosus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Desulfurococcus_mucosus], "data/meta_species/species/Desulfurococcus_mucosus.fa");
  panel_file_paths[Dialister_invisus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Dialister_invisus], "data/meta_species/species/Dialister_invisus.fa");
  panel_file_paths[Dialister_microaerophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Dialister_microaerophilus], "data/meta_species/species/Dialister_microaerophilus.fa");
  panel_file_paths[Dickeya_dadantii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Dickeya_dadantii], "data/meta_species/species/Dickeya_dadantii.fa");
  panel_file_paths[Dickeya_zeae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Dickeya_zeae], "data/meta_species/species/Dickeya_zeae.fa");
  panel_file_paths[Dictyoglomus_thermophilum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Dictyoglomus_thermophilum], "data/meta_species/species/Dictyoglomus_thermophilum.fa");
  panel_file_paths[Dictyoglomus_turgidum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Dictyoglomus_turgidum], "data/meta_species/species/Dictyoglomus_turgidum.fa");
  panel_file_paths[Dorea_formicigenerans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Dorea_formicigenerans], "data/meta_species/species/Dorea_formicigenerans.fa");
  panel_file_paths[Dorea_longicatena] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Dorea_longicatena], "data/meta_species/species/Dorea_longicatena.fa");
  panel_file_paths[Edwardsiella_ictaluri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Edwardsiella_ictaluri], "data/meta_species/species/Edwardsiella_ictaluri.fa");
  panel_file_paths[Edwardsiella_tarda] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Edwardsiella_tarda], "data/meta_species/species/Edwardsiella_tarda.fa");
  panel_file_paths[Eggerthella_lenta] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eggerthella_lenta], "data/meta_species/species/Eggerthella_lenta.fa");
  panel_file_paths[Ehrlichia_canis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ehrlichia_canis], "data/meta_species/species/Ehrlichia_canis.fa");
  panel_file_paths[Ehrlichia_chaffeensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ehrlichia_chaffeensis], "data/meta_species/species/Ehrlichia_chaffeensis.fa");
  panel_file_paths[Ehrlichia_ruminantium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ehrlichia_ruminantium], "data/meta_species/species/Ehrlichia_ruminantium.fa");
  panel_file_paths[Ensifer_medicae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ensifer_medicae], "data/meta_species/species/Ensifer_medicae.fa");
  panel_file_paths[Ensifer_meliloti] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ensifer_meliloti], "data/meta_species/species/Ensifer_meliloti.fa");
  panel_file_paths[Enterobacter_cancerogenus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Enterobacter_cancerogenus], "data/meta_species/species/Enterobacter_cancerogenus.fa");
  panel_file_paths[Enterobacter_cloacae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Enterobacter_cloacae], "data/meta_species/species/Enterobacter_cloacae.fa");
  panel_file_paths[Enterococcus_casseliflavus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Enterococcus_casseliflavus], "data/meta_species/species/Enterococcus_casseliflavus.fa");
  panel_file_paths[Enterococcus_faecalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Enterococcus_faecalis], "data/meta_species/species/Enterococcus_faecalis.fa");
  panel_file_paths[Enterococcus_faecium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Enterococcus_faecium], "data/meta_species/species/Enterococcus_faecium.fa");
  panel_file_paths[Enterococcus_gallinarum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Enterococcus_gallinarum], "data/meta_species/species/Enterococcus_gallinarum.fa");
  panel_file_paths[Enterococcus_italicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Enterococcus_italicus], "data/meta_species/species/Enterococcus_italicus.fa");
  panel_file_paths[Erwinia_amylovora] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Erwinia_amylovora], "data/meta_species/species/Erwinia_amylovora.fa");
  panel_file_paths[Erwinia_billingiae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Erwinia_billingiae], "data/meta_species/species/Erwinia_billingiae.fa");
  panel_file_paths[Erwinia_pyrifoliae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Erwinia_pyrifoliae], "data/meta_species/species/Erwinia_pyrifoliae.fa");
  panel_file_paths[Erwinia_tasmaniensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Erwinia_tasmaniensis], "data/meta_species/species/Erwinia_tasmaniensis.fa");
  panel_file_paths[Erythrobacter_litoralis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Erythrobacter_litoralis], "data/meta_species/species/Erythrobacter_litoralis.fa");
  panel_file_paths[Escherichia_albertii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Escherichia_albertii], "data/meta_species/species/Escherichia_albertii.fa");
  panel_file_paths[Escherichia_coli] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Escherichia_coli], "data/meta_species/species/Escherichia_coli.fa");
  panel_file_paths[Escherichia_fergusonii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Escherichia_fergusonii], "data/meta_species/species/Escherichia_fergusonii.fa");
  panel_file_paths[Eubacterium_cellulosolvens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eubacterium_cellulosolvens], "data/meta_species/species/Eubacterium_cellulosolvens.fa");
  panel_file_paths[Eubacterium_eligens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eubacterium_eligens], "data/meta_species/species/Eubacterium_eligens.fa");
  panel_file_paths[Eubacterium_hallii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eubacterium_hallii], "data/meta_species/species/Eubacterium_hallii.fa");
  panel_file_paths[Eubacterium_limosum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eubacterium_limosum], "data/meta_species/species/Eubacterium_limosum.fa");
  panel_file_paths[Eubacterium_rectale] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eubacterium_rectale], "data/meta_species/species/Eubacterium_rectale.fa");
  panel_file_paths[Eubacterium_saburreum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eubacterium_saburreum], "data/meta_species/species/Eubacterium_saburreum.fa");
  panel_file_paths[Eubacterium_saphenum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eubacterium_saphenum], "data/meta_species/species/Eubacterium_saphenum.fa");
  panel_file_paths[Eubacterium_siraeum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eubacterium_siraeum], "data/meta_species/species/Eubacterium_siraeum.fa");
  panel_file_paths[Eubacterium_ventriosum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Eubacterium_ventriosum], "data/meta_species/species/Eubacterium_ventriosum.fa");
  panel_file_paths[Exiguobacterium_sibiricum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Exiguobacterium_sibiricum], "data/meta_species/species/Exiguobacterium_sibiricum.fa");
  panel_file_paths[Faecalibacterium_cf] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Faecalibacterium_cf], "data/meta_species/species/Faecalibacterium_cf.fa");
  panel_file_paths[Faecalibacterium_prausnitzii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Faecalibacterium_prausnitzii], "data/meta_species/species/Faecalibacterium_prausnitzii.fa");
  panel_file_paths[Flavobacterium_johnsoniae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Flavobacterium_johnsoniae], "data/meta_species/species/Flavobacterium_johnsoniae.fa");
  panel_file_paths[Flavobacterium_psychrophilum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Flavobacterium_psychrophilum], "data/meta_species/species/Flavobacterium_psychrophilum.fa");
  panel_file_paths[Francisella_novicida] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Francisella_novicida], "data/meta_species/species/Francisella_novicida.fa");
  panel_file_paths[Francisella_philomiragia] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Francisella_philomiragia], "data/meta_species/species/Francisella_philomiragia.fa");
  panel_file_paths[Francisella_tularensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Francisella_tularensis], "data/meta_species/species/Francisella_tularensis.fa");
  panel_file_paths[Frankia_alni] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Frankia_alni], "data/meta_species/species/Frankia_alni.fa");
  panel_file_paths[Frankia_symbiont] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Frankia_symbiont], "data/meta_species/species/Frankia_symbiont.fa");
  panel_file_paths[Fusobacterium_gonidiaformans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Fusobacterium_gonidiaformans], "data/meta_species/species/Fusobacterium_gonidiaformans.fa");
  panel_file_paths[Fusobacterium_mortiferum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Fusobacterium_mortiferum], "data/meta_species/species/Fusobacterium_mortiferum.fa");
  panel_file_paths[Fusobacterium_nucleatum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Fusobacterium_nucleatum], "data/meta_species/species/Fusobacterium_nucleatum.fa");
  panel_file_paths[Fusobacterium_periodonticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Fusobacterium_periodonticum], "data/meta_species/species/Fusobacterium_periodonticum.fa");
  panel_file_paths[Fusobacterium_ulcerans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Fusobacterium_ulcerans], "data/meta_species/species/Fusobacterium_ulcerans.fa");
  panel_file_paths[Fusobacterium_varium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Fusobacterium_varium], "data/meta_species/species/Fusobacterium_varium.fa");
  panel_file_paths[Gemella_haemolysans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Gemella_haemolysans], "data/meta_species/species/Gemella_haemolysans.fa");
  panel_file_paths[Gemella_moribillum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Gemella_moribillum], "data/meta_species/species/Gemella_moribillum.fa");
  panel_file_paths[Geobacillus_kaustophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Geobacillus_kaustophilus], "data/meta_species/species/Geobacillus_kaustophilus.fa");
  panel_file_paths[Geobacillus_thermodenitrificans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Geobacillus_thermodenitrificans], "data/meta_species/species/Geobacillus_thermodenitrificans.fa");
  panel_file_paths[Geobacillus_thermoglucosidasius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Geobacillus_thermoglucosidasius], "data/meta_species/species/Geobacillus_thermoglucosidasius.fa");
  panel_file_paths[Geobacter_bemidjiensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Geobacter_bemidjiensis], "data/meta_species/species/Geobacter_bemidjiensis.fa");
  panel_file_paths[Geobacter_lovleyi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Geobacter_lovleyi], "data/meta_species/species/Geobacter_lovleyi.fa");
  panel_file_paths[Geobacter_metallireducens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Geobacter_metallireducens], "data/meta_species/species/Geobacter_metallireducens.fa");
  panel_file_paths[Geobacter_sulfurreducens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Geobacter_sulfurreducens], "data/meta_species/species/Geobacter_sulfurreducens.fa");
  panel_file_paths[Geobacter_uraniumreducens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Geobacter_uraniumreducens], "data/meta_species/species/Geobacter_uraniumreducens.fa");
  panel_file_paths[Gluconacetobacter_diazotrophicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Gluconacetobacter_diazotrophicus], "data/meta_species/species/Gluconacetobacter_diazotrophicus.fa");
  panel_file_paths[Gluconacetobacter_hansenii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Gluconacetobacter_hansenii], "data/meta_species/species/Gluconacetobacter_hansenii.fa");
  panel_file_paths[Granulicatella_adiacens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Granulicatella_adiacens], "data/meta_species/species/Granulicatella_adiacens.fa");
  panel_file_paths[Granulicatella_elegans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Granulicatella_elegans], "data/meta_species/species/Granulicatella_elegans.fa");
  panel_file_paths[Haemophilus_ducreyi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Haemophilus_ducreyi], "data/meta_species/species/Haemophilus_ducreyi.fa");
  panel_file_paths[Haemophilus_influenzae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Haemophilus_influenzae], "data/meta_species/species/Haemophilus_influenzae.fa");
  panel_file_paths[Haemophilus_parainfluenzae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Haemophilus_parainfluenzae], "data/meta_species/species/Haemophilus_parainfluenzae.fa");
  panel_file_paths[Haemophilus_parasuis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Haemophilus_parasuis], "data/meta_species/species/Haemophilus_parasuis.fa");
  panel_file_paths[Halanaerobium_praevalens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Halanaerobium_praevalens], "data/meta_species/species/Halanaerobium_praevalens.fa");
  panel_file_paths[Halobacterium_salinarum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Halobacterium_salinarum], "data/meta_species/species/Halobacterium_salinarum.fa");
  panel_file_paths[Helicobacter_acinonychis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_acinonychis], "data/meta_species/species/Helicobacter_acinonychis.fa");
  panel_file_paths[Helicobacter_bilis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_bilis], "data/meta_species/species/Helicobacter_bilis.fa");
  panel_file_paths[Helicobacter_canadensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_canadensis], "data/meta_species/species/Helicobacter_canadensis.fa");
  panel_file_paths[Helicobacter_cinaedi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_cinaedi], "data/meta_species/species/Helicobacter_cinaedi.fa");
  panel_file_paths[Helicobacter_felis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_felis], "data/meta_species/species/Helicobacter_felis.fa");
  panel_file_paths[Helicobacter_hepaticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_hepaticus], "data/meta_species/species/Helicobacter_hepaticus.fa");
  panel_file_paths[Helicobacter_mustelae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_mustelae], "data/meta_species/species/Helicobacter_mustelae.fa");
  panel_file_paths[Helicobacter_pullorum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_pullorum], "data/meta_species/species/Helicobacter_pullorum.fa");
  panel_file_paths[Helicobacter_pylori] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_pylori], "data/meta_species/species/Helicobacter_pylori.fa");
  panel_file_paths[Helicobacter_suis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_suis], "data/meta_species/species/Helicobacter_suis.fa");
  panel_file_paths[Helicobacter_winghamensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Helicobacter_winghamensis], "data/meta_species/species/Helicobacter_winghamensis.fa");
  panel_file_paths[Idiomarina_baltica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Idiomarina_baltica], "data/meta_species/species/Idiomarina_baltica.fa");
  panel_file_paths[Idiomarina_loihiensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Idiomarina_loihiensis], "data/meta_species/species/Idiomarina_loihiensis.fa");
  panel_file_paths[Kingella_denitrificans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Kingella_denitrificans], "data/meta_species/species/Kingella_denitrificans.fa");
  panel_file_paths[Kingella_oralis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Kingella_oralis], "data/meta_species/species/Kingella_oralis.fa");
  panel_file_paths[Klebsiella_pneumoniae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Klebsiella_pneumoniae], "data/meta_species/species/Klebsiella_pneumoniae.fa");
  panel_file_paths[Klebsiella_variicola] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Klebsiella_variicola], "data/meta_species/species/Klebsiella_variicola.fa");
  panel_file_paths[Labrenzia_aggregata] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Labrenzia_aggregata], "data/meta_species/species/Labrenzia_aggregata.fa");
  panel_file_paths[Labrenzia_alexandrii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Labrenzia_alexandrii], "data/meta_species/species/Labrenzia_alexandrii.fa");
  panel_file_paths[Lactobacillus_acidophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_acidophilus], "data/meta_species/species/Lactobacillus_acidophilus.fa");
  panel_file_paths[Lactobacillus_amylolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_amylolyticus], "data/meta_species/species/Lactobacillus_amylolyticus.fa");
  panel_file_paths[Lactobacillus_amylovorus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_amylovorus], "data/meta_species/species/Lactobacillus_amylovorus.fa");
  panel_file_paths[Lactobacillus_antri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_antri], "data/meta_species/species/Lactobacillus_antri.fa");
  panel_file_paths[Lactobacillus_brevis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_brevis], "data/meta_species/species/Lactobacillus_brevis.fa");
  panel_file_paths[Lactobacillus_buchneri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_buchneri], "data/meta_species/species/Lactobacillus_buchneri.fa");
  panel_file_paths[Lactobacillus_casei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_casei], "data/meta_species/species/Lactobacillus_casei.fa");
  panel_file_paths[Lactobacillus_coleohominis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_coleohominis], "data/meta_species/species/Lactobacillus_coleohominis.fa");
  panel_file_paths[Lactobacillus_crispatus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_crispatus], "data/meta_species/species/Lactobacillus_crispatus.fa");
  panel_file_paths[Lactobacillus_delbrueckii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_delbrueckii], "data/meta_species/species/Lactobacillus_delbrueckii.fa");
  panel_file_paths[Lactobacillus_fermentum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_fermentum], "data/meta_species/species/Lactobacillus_fermentum.fa");
  panel_file_paths[Lactobacillus_gasseri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_gasseri], "data/meta_species/species/Lactobacillus_gasseri.fa");
  panel_file_paths[Lactobacillus_helveticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_helveticus], "data/meta_species/species/Lactobacillus_helveticus.fa");
  panel_file_paths[Lactobacillus_hilgardii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_hilgardii], "data/meta_species/species/Lactobacillus_hilgardii.fa");
  panel_file_paths[Lactobacillus_iners] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_iners], "data/meta_species/species/Lactobacillus_iners.fa");
  panel_file_paths[Lactobacillus_jensenii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_jensenii], "data/meta_species/species/Lactobacillus_jensenii.fa");
  panel_file_paths[Lactobacillus_johnsonii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_johnsonii], "data/meta_species/species/Lactobacillus_johnsonii.fa");
  panel_file_paths[Lactobacillus_oris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_oris], "data/meta_species/species/Lactobacillus_oris.fa");
  panel_file_paths[Lactobacillus_paracasei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_paracasei], "data/meta_species/species/Lactobacillus_paracasei.fa");
  panel_file_paths[Lactobacillus_plantarum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_plantarum], "data/meta_species/species/Lactobacillus_plantarum.fa");
  panel_file_paths[Lactobacillus_reuteri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_reuteri], "data/meta_species/species/Lactobacillus_reuteri.fa");
  panel_file_paths[Lactobacillus_rhamnosus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_rhamnosus], "data/meta_species/species/Lactobacillus_rhamnosus.fa");
  panel_file_paths[Lactobacillus_ruminis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_ruminis], "data/meta_species/species/Lactobacillus_ruminis.fa");
  panel_file_paths[Lactobacillus_sakei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_sakei], "data/meta_species/species/Lactobacillus_sakei.fa");
  panel_file_paths[Lactobacillus_salivarius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_salivarius], "data/meta_species/species/Lactobacillus_salivarius.fa");
  panel_file_paths[Lactobacillus_ultunensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_ultunensis], "data/meta_species/species/Lactobacillus_ultunensis.fa");
  panel_file_paths[Lactobacillus_vaginalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lactobacillus_vaginalis], "data/meta_species/species/Lactobacillus_vaginalis.fa");
  panel_file_paths[Legionella_drancourtii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Legionella_drancourtii], "data/meta_species/species/Legionella_drancourtii.fa");
  panel_file_paths[Legionella_longbeachae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Legionella_longbeachae], "data/meta_species/species/Legionella_longbeachae.fa");
  panel_file_paths[Legionella_pneumophila] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Legionella_pneumophila], "data/meta_species/species/Legionella_pneumophila.fa");
  panel_file_paths[Leptospira_biflexa] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leptospira_biflexa], "data/meta_species/species/Leptospira_biflexa.fa");
  panel_file_paths[Leptospira_borgpetersenii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leptospira_borgpetersenii], "data/meta_species/species/Leptospira_borgpetersenii.fa");
  panel_file_paths[Leptospira_interrogans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leptospira_interrogans], "data/meta_species/species/Leptospira_interrogans.fa");
  panel_file_paths[Leptotrichia_buccalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leptotrichia_buccalis], "data/meta_species/species/Leptotrichia_buccalis.fa");
  panel_file_paths[Leptotrichia_goodfellowii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leptotrichia_goodfellowii], "data/meta_species/species/Leptotrichia_goodfellowii.fa");
  panel_file_paths[Leptotrichia_hofstadii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leptotrichia_hofstadii], "data/meta_species/species/Leptotrichia_hofstadii.fa");
  panel_file_paths[Leuconostoc_citreum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leuconostoc_citreum], "data/meta_species/species/Leuconostoc_citreum.fa");
  panel_file_paths[Leuconostoc_gasicomitatum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leuconostoc_gasicomitatum], "data/meta_species/species/Leuconostoc_gasicomitatum.fa");
  panel_file_paths[Leuconostoc_kimchii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leuconostoc_kimchii], "data/meta_species/species/Leuconostoc_kimchii.fa");
  panel_file_paths[Leuconostoc_mesenteroides] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Leuconostoc_mesenteroides], "data/meta_species/species/Leuconostoc_mesenteroides.fa");
  panel_file_paths[Listeria_grayi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Listeria_grayi], "data/meta_species/species/Listeria_grayi.fa");
  panel_file_paths[Listeria_innocua] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Listeria_innocua], "data/meta_species/species/Listeria_innocua.fa");
  panel_file_paths[Listeria_monocytogenes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Listeria_monocytogenes], "data/meta_species/species/Listeria_monocytogenes.fa");
  panel_file_paths[Listeria_seeligeri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Listeria_seeligeri], "data/meta_species/species/Listeria_seeligeri.fa");
  panel_file_paths[Listeria_welshimeri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Listeria_welshimeri], "data/meta_species/species/Listeria_welshimeri.fa");
  panel_file_paths[Loktanella_vestfoldensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Loktanella_vestfoldensis], "data/meta_species/species/Loktanella_vestfoldensis.fa");
  panel_file_paths[Lysinibacillus_fusiformis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lysinibacillus_fusiformis], "data/meta_species/species/Lysinibacillus_fusiformis.fa");
  panel_file_paths[Lysinibacillus_sphaericus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Lysinibacillus_sphaericus], "data/meta_species/species/Lysinibacillus_sphaericus.fa");
  panel_file_paths[Magnetospirillum_magneticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Magnetospirillum_magneticum], "data/meta_species/species/Magnetospirillum_magneticum.fa");
  panel_file_paths[Magnetospirillum_magnetotacticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Magnetospirillum_magnetotacticum], "data/meta_species/species/Magnetospirillum_magnetotacticum.fa");
  panel_file_paths[Marinobacter_algicola] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Marinobacter_algicola], "data/meta_species/species/Marinobacter_algicola.fa");
  panel_file_paths[Marinobacter_aquaeolei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Marinobacter_aquaeolei], "data/meta_species/species/Marinobacter_aquaeolei.fa");
  panel_file_paths[Marinobacter_bacterium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Marinobacter_bacterium], "data/meta_species/species/Marinobacter_bacterium.fa");
  panel_file_paths[Megasphaera_micronuciformis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Megasphaera_micronuciformis], "data/meta_species/species/Megasphaera_micronuciformis.fa");
  panel_file_paths[Meiothermus_ruber] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Meiothermus_ruber], "data/meta_species/species/Meiothermus_ruber.fa");
  panel_file_paths[Meiothermus_silvanus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Meiothermus_silvanus], "data/meta_species/species/Meiothermus_silvanus.fa");
  panel_file_paths[Mesorhizobium_ciceri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mesorhizobium_ciceri], "data/meta_species/species/Mesorhizobium_ciceri.fa");
  panel_file_paths[Mesorhizobium_loti] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mesorhizobium_loti], "data/meta_species/species/Mesorhizobium_loti.fa");
  panel_file_paths[Mesorhizobium_opportunistum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mesorhizobium_opportunistum], "data/meta_species/species/Mesorhizobium_opportunistum.fa");
  panel_file_paths[Methanobrevibacter_ruminantium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanobrevibacter_ruminantium], "data/meta_species/species/Methanobrevibacter_ruminantium.fa");
  panel_file_paths[Methanobrevibacter_smithii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanobrevibacter_smithii], "data/meta_species/species/Methanobrevibacter_smithii.fa");
  panel_file_paths[Methanocaldococcus_fervens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanocaldococcus_fervens], "data/meta_species/species/Methanocaldococcus_fervens.fa");
  panel_file_paths[Methanocaldococcus_infernus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanocaldococcus_infernus], "data/meta_species/species/Methanocaldococcus_infernus.fa");
  panel_file_paths[Methanocaldococcus_jannaschii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanocaldococcus_jannaschii], "data/meta_species/species/Methanocaldococcus_jannaschii.fa");
  panel_file_paths[Methanocaldococcus_vulcanius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanocaldococcus_vulcanius], "data/meta_species/species/Methanocaldococcus_vulcanius.fa");
  panel_file_paths[Methanocella_paludicola] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanocella_paludicola], "data/meta_species/species/Methanocella_paludicola.fa");
  panel_file_paths[Methanococcus_aeolicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanococcus_aeolicus], "data/meta_species/species/Methanococcus_aeolicus.fa");
  panel_file_paths[Methanococcus_maripaludis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanococcus_maripaludis], "data/meta_species/species/Methanococcus_maripaludis.fa");
  panel_file_paths[Methanococcus_vannielii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanococcus_vannielii], "data/meta_species/species/Methanococcus_vannielii.fa");
  panel_file_paths[Methanococcus_voltae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanococcus_voltae], "data/meta_species/species/Methanococcus_voltae.fa");
  panel_file_paths[Methanosarcina_acetivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanosarcina_acetivorans], "data/meta_species/species/Methanosarcina_acetivorans.fa");
  panel_file_paths[Methanosarcina_barkeri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanosarcina_barkeri], "data/meta_species/species/Methanosarcina_barkeri.fa");
  panel_file_paths[Methanosarcina_mazei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanosarcina_mazei], "data/meta_species/species/Methanosarcina_mazei.fa");
  panel_file_paths[Methanothermobacter_marburgensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanothermobacter_marburgensis], "data/meta_species/species/Methanothermobacter_marburgensis.fa");
  panel_file_paths[Methanothermobacter_thermautotrophicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methanothermobacter_thermautotrophicus], "data/meta_species/species/Methanothermobacter_thermautotrophicus.fa");
  panel_file_paths[Methylobacterium_chloromethanicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methylobacterium_chloromethanicum], "data/meta_species/species/Methylobacterium_chloromethanicum.fa");
  panel_file_paths[Methylobacterium_extorquens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methylobacterium_extorquens], "data/meta_species/species/Methylobacterium_extorquens.fa");
  panel_file_paths[Methylobacterium_nodulans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methylobacterium_nodulans], "data/meta_species/species/Methylobacterium_nodulans.fa");
  panel_file_paths[Methylobacterium_populi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methylobacterium_populi], "data/meta_species/species/Methylobacterium_populi.fa");
  panel_file_paths[Methylobacterium_radiotolerans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methylobacterium_radiotolerans], "data/meta_species/species/Methylobacterium_radiotolerans.fa");
  panel_file_paths[Methylotenera_mobilis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Methylotenera_mobilis], "data/meta_species/species/Methylotenera_mobilis.fa");
  panel_file_paths[Micromonospora_aurantiaca] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Micromonospora_aurantiaca], "data/meta_species/species/Micromonospora_aurantiaca.fa");
  panel_file_paths[Mobiluncus_curtisii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mobiluncus_curtisii], "data/meta_species/species/Mobiluncus_curtisii.fa");
  panel_file_paths[Mobiluncus_mulieris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mobiluncus_mulieris], "data/meta_species/species/Mobiluncus_mulieris.fa");
  panel_file_paths[Mycobacterium_abscessus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_abscessus], "data/meta_species/species/Mycobacterium_abscessus.fa");
  panel_file_paths[Mycobacterium_avium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_avium], "data/meta_species/species/Mycobacterium_avium.fa");
  panel_file_paths[Mycobacterium_bovis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_bovis], "data/meta_species/species/Mycobacterium_bovis.fa");
  panel_file_paths[Mycobacterium_gilvum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_gilvum], "data/meta_species/species/Mycobacterium_gilvum.fa");
  panel_file_paths[Mycobacterium_intracellulare] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_intracellulare], "data/meta_species/species/Mycobacterium_intracellulare.fa");
  panel_file_paths[Mycobacterium_kansasii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_kansasii], "data/meta_species/species/Mycobacterium_kansasii.fa");
  panel_file_paths[Mycobacterium_leprae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_leprae], "data/meta_species/species/Mycobacterium_leprae.fa");
  panel_file_paths[Mycobacterium_marinum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_marinum], "data/meta_species/species/Mycobacterium_marinum.fa");
  panel_file_paths[Mycobacterium_parascrofulaceum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_parascrofulaceum], "data/meta_species/species/Mycobacterium_parascrofulaceum.fa");
  panel_file_paths[Mycobacterium_smegmatis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_smegmatis], "data/meta_species/species/Mycobacterium_smegmatis.fa");
  panel_file_paths[Mycobacterium_tuberculosis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_tuberculosis], "data/meta_species/species/Mycobacterium_tuberculosis.fa");
  panel_file_paths[Mycobacterium_ulcerans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_ulcerans], "data/meta_species/species/Mycobacterium_ulcerans.fa");
  panel_file_paths[Mycobacterium_vanbaalenii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycobacterium_vanbaalenii], "data/meta_species/species/Mycobacterium_vanbaalenii.fa");
  panel_file_paths[Mycoplasma_agalactiae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_agalactiae], "data/meta_species/species/Mycoplasma_agalactiae.fa");
  panel_file_paths[Mycoplasma_alligatoris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_alligatoris], "data/meta_species/species/Mycoplasma_alligatoris.fa");
  panel_file_paths[Mycoplasma_arthritidis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_arthritidis], "data/meta_species/species/Mycoplasma_arthritidis.fa");
  panel_file_paths[Mycoplasma_bovis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_bovis], "data/meta_species/species/Mycoplasma_bovis.fa");
  panel_file_paths[Mycoplasma_capricolum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_capricolum], "data/meta_species/species/Mycoplasma_capricolum.fa");
  panel_file_paths[Mycoplasma_conjunctivae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_conjunctivae], "data/meta_species/species/Mycoplasma_conjunctivae.fa");
  panel_file_paths[Mycoplasma_crocodyli] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_crocodyli], "data/meta_species/species/Mycoplasma_crocodyli.fa");
  panel_file_paths[Mycoplasma_fermentans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_fermentans], "data/meta_species/species/Mycoplasma_fermentans.fa");
  panel_file_paths[Mycoplasma_gallisepticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_gallisepticum], "data/meta_species/species/Mycoplasma_gallisepticum.fa");
  panel_file_paths[Mycoplasma_genitalium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_genitalium], "data/meta_species/species/Mycoplasma_genitalium.fa");
  panel_file_paths[Mycoplasma_haemofelis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_haemofelis], "data/meta_species/species/Mycoplasma_haemofelis.fa");
  panel_file_paths[Mycoplasma_hominis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_hominis], "data/meta_species/species/Mycoplasma_hominis.fa");
  panel_file_paths[Mycoplasma_hyopneumoniae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_hyopneumoniae], "data/meta_species/species/Mycoplasma_hyopneumoniae.fa");
  panel_file_paths[Mycoplasma_hyorhinis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_hyorhinis], "data/meta_species/species/Mycoplasma_hyorhinis.fa");
  panel_file_paths[Mycoplasma_leachii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_leachii], "data/meta_species/species/Mycoplasma_leachii.fa");
  panel_file_paths[Mycoplasma_mobile] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_mobile], "data/meta_species/species/Mycoplasma_mobile.fa");
  panel_file_paths[Mycoplasma_mycoides] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_mycoides], "data/meta_species/species/Mycoplasma_mycoides.fa");
  panel_file_paths[Mycoplasma_penetrans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_penetrans], "data/meta_species/species/Mycoplasma_penetrans.fa");
  panel_file_paths[Mycoplasma_pneumoniae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_pneumoniae], "data/meta_species/species/Mycoplasma_pneumoniae.fa");
  panel_file_paths[Mycoplasma_pulmonis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_pulmonis], "data/meta_species/species/Mycoplasma_pulmonis.fa");
  panel_file_paths[Mycoplasma_synoviae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Mycoplasma_synoviae], "data/meta_species/species/Mycoplasma_synoviae.fa");
  panel_file_paths[Neisseria_cinerea] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_cinerea], "data/meta_species/species/Neisseria_cinerea.fa");
  panel_file_paths[Neisseria_elongata] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_elongata], "data/meta_species/species/Neisseria_elongata.fa");
  panel_file_paths[Neisseria_flavescens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_flavescens], "data/meta_species/species/Neisseria_flavescens.fa");
  panel_file_paths[Neisseria_gonorrhoeae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_gonorrhoeae], "data/meta_species/species/Neisseria_gonorrhoeae.fa");
  panel_file_paths[Neisseria_lactamica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_lactamica], "data/meta_species/species/Neisseria_lactamica.fa");
  panel_file_paths[Neisseria_meningitidis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_meningitidis], "data/meta_species/species/Neisseria_meningitidis.fa");
  panel_file_paths[Neisseria_mucosa] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_mucosa], "data/meta_species/species/Neisseria_mucosa.fa");
  panel_file_paths[Neisseria_polysaccharea] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_polysaccharea], "data/meta_species/species/Neisseria_polysaccharea.fa");
  panel_file_paths[Neisseria_sicca] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_sicca], "data/meta_species/species/Neisseria_sicca.fa");
  panel_file_paths[Neisseria_subflava] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neisseria_subflava], "data/meta_species/species/Neisseria_subflava.fa");
  panel_file_paths[Neorickettsia_risticii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neorickettsia_risticii], "data/meta_species/species/Neorickettsia_risticii.fa");
  panel_file_paths[Neorickettsia_sennetsu] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Neorickettsia_sennetsu], "data/meta_species/species/Neorickettsia_sennetsu.fa");
  panel_file_paths[Nitrobacter_hamburgensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Nitrobacter_hamburgensis], "data/meta_species/species/Nitrobacter_hamburgensis.fa");
  panel_file_paths[Nitrobacter_winogradskyi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Nitrobacter_winogradskyi], "data/meta_species/species/Nitrobacter_winogradskyi.fa");
  panel_file_paths[Nitrosococcus_halophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Nitrosococcus_halophilus], "data/meta_species/species/Nitrosococcus_halophilus.fa");
  panel_file_paths[Nitrosococcus_oceani] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Nitrosococcus_oceani], "data/meta_species/species/Nitrosococcus_oceani.fa");
  panel_file_paths[Nitrosococcus_watsoni] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Nitrosococcus_watsoni], "data/meta_species/species/Nitrosococcus_watsoni.fa");
  panel_file_paths[Nitrosomonas_europaea] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Nitrosomonas_europaea], "data/meta_species/species/Nitrosomonas_europaea.fa");
  panel_file_paths[Nitrosomonas_eutropha] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Nitrosomonas_eutropha], "data/meta_species/species/Nitrosomonas_eutropha.fa");
  panel_file_paths[Nostoc_punctiforme] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Nostoc_punctiforme], "data/meta_species/species/Nostoc_punctiforme.fa");
  panel_file_paths[Oceanicola_batsensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Oceanicola_batsensis], "data/meta_species/species/Oceanicola_batsensis.fa");
  panel_file_paths[Oceanicola_granulosus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Oceanicola_granulosus], "data/meta_species/species/Oceanicola_granulosus.fa");
  panel_file_paths[Ochrobactrum_anthropi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ochrobactrum_anthropi], "data/meta_species/species/Ochrobactrum_anthropi.fa");
  panel_file_paths[Ochrobactrum_intermedium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ochrobactrum_intermedium], "data/meta_species/species/Ochrobactrum_intermedium.fa");
  panel_file_paths[Oribacterium_sinus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Oribacterium_sinus], "data/meta_species/species/Oribacterium_sinus.fa");
  panel_file_paths[Paenibacillus_curdlanolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Paenibacillus_curdlanolyticus], "data/meta_species/species/Paenibacillus_curdlanolyticus.fa");
  panel_file_paths[Paenibacillus_larvae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Paenibacillus_larvae], "data/meta_species/species/Paenibacillus_larvae.fa");
  panel_file_paths[Paenibacillus_polymyxa] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Paenibacillus_polymyxa], "data/meta_species/species/Paenibacillus_polymyxa.fa");
  panel_file_paths[Paenibacillus_vortex] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Paenibacillus_vortex], "data/meta_species/species/Paenibacillus_vortex.fa");
  panel_file_paths[Pantoea_ananatis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pantoea_ananatis], "data/meta_species/species/Pantoea_ananatis.fa");
  panel_file_paths[Pantoea_vagans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pantoea_vagans], "data/meta_species/species/Pantoea_vagans.fa");
  panel_file_paths[Parabacteroides_distasonis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Parabacteroides_distasonis], "data/meta_species/species/Parabacteroides_distasonis.fa");
  panel_file_paths[Parabacteroides_johnsonii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Parabacteroides_johnsonii], "data/meta_species/species/Parabacteroides_johnsonii.fa");
  panel_file_paths[Parabacteroides_merdae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Parabacteroides_merdae], "data/meta_species/species/Parabacteroides_merdae.fa");
  panel_file_paths[Pasteurella_dagmatis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pasteurella_dagmatis], "data/meta_species/species/Pasteurella_dagmatis.fa");
  panel_file_paths[Pasteurella_multocida] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pasteurella_multocida], "data/meta_species/species/Pasteurella_multocida.fa");
  panel_file_paths[Pectobacterium_carotovora] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pectobacterium_carotovora], "data/meta_species/species/Pectobacterium_carotovora.fa");
  panel_file_paths[Pectobacterium_carotovorum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pectobacterium_carotovorum], "data/meta_species/species/Pectobacterium_carotovorum.fa");
  panel_file_paths[Pectobacterium_wasabiae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pectobacterium_wasabiae], "data/meta_species/species/Pectobacterium_wasabiae.fa");
  panel_file_paths[Pediococcus_acidilactici] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pediococcus_acidilactici], "data/meta_species/species/Pediococcus_acidilactici.fa");
  panel_file_paths[Pediococcus_pentosaceus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pediococcus_pentosaceus], "data/meta_species/species/Pediococcus_pentosaceus.fa");
  panel_file_paths[Pedobacter_heparinus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pedobacter_heparinus], "data/meta_species/species/Pedobacter_heparinus.fa");
  panel_file_paths[Pedobacter_saltans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pedobacter_saltans], "data/meta_species/species/Pedobacter_saltans.fa");
  panel_file_paths[Pelobacter_carbinolicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pelobacter_carbinolicus], "data/meta_species/species/Pelobacter_carbinolicus.fa");
  panel_file_paths[Pelobacter_propionicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pelobacter_propionicus], "data/meta_species/species/Pelobacter_propionicus.fa");
  panel_file_paths[Pelodictyon_clathratiforme] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pelodictyon_clathratiforme], "data/meta_species/species/Pelodictyon_clathratiforme.fa");
  panel_file_paths[Pelodictyon_luteolum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pelodictyon_luteolum], "data/meta_species/species/Pelodictyon_luteolum.fa");
  panel_file_paths[Peptoniphilus_duerdenii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Peptoniphilus_duerdenii], "data/meta_species/species/Peptoniphilus_duerdenii.fa");
  panel_file_paths[Peptoniphilus_harei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Peptoniphilus_harei], "data/meta_species/species/Peptoniphilus_harei.fa");
  panel_file_paths[Peptoniphilus_lacrimalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Peptoniphilus_lacrimalis], "data/meta_species/species/Peptoniphilus_lacrimalis.fa");
  panel_file_paths[Peptostreptococcus_anaerobius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Peptostreptococcus_anaerobius], "data/meta_species/species/Peptostreptococcus_anaerobius.fa");
  panel_file_paths[Peptostreptococcus_stomatis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Peptostreptococcus_stomatis], "data/meta_species/species/Peptostreptococcus_stomatis.fa");
  panel_file_paths[Photobacterium_angustum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Photobacterium_angustum], "data/meta_species/species/Photobacterium_angustum.fa");
  panel_file_paths[Photobacterium_damselae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Photobacterium_damselae], "data/meta_species/species/Photobacterium_damselae.fa");
  panel_file_paths[Photobacterium_profundum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Photobacterium_profundum], "data/meta_species/species/Photobacterium_profundum.fa");
  panel_file_paths[Photorhabdus_asymbiotica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Photorhabdus_asymbiotica], "data/meta_species/species/Photorhabdus_asymbiotica.fa");
  panel_file_paths[Photorhabdus_luminescens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Photorhabdus_luminescens], "data/meta_species/species/Photorhabdus_luminescens.fa");
  panel_file_paths[Planctomyces_brasiliensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Planctomyces_brasiliensis], "data/meta_species/species/Planctomyces_brasiliensis.fa");
  panel_file_paths[Planctomyces_limnophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Planctomyces_limnophilus], "data/meta_species/species/Planctomyces_limnophilus.fa");
  panel_file_paths[Planctomyces_maris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Planctomyces_maris], "data/meta_species/species/Planctomyces_maris.fa");
  panel_file_paths[Polaribacter_irgensii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Polaribacter_irgensii], "data/meta_species/species/Polaribacter_irgensii.fa");
  panel_file_paths[Polaromonas_naphthalenivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Polaromonas_naphthalenivorans], "data/meta_species/species/Polaromonas_naphthalenivorans.fa");
  panel_file_paths[Polynucleobacter_necessarius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Polynucleobacter_necessarius], "data/meta_species/species/Polynucleobacter_necessarius.fa");
  panel_file_paths[Porphyromonas_asaccharolytica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Porphyromonas_asaccharolytica], "data/meta_species/species/Porphyromonas_asaccharolytica.fa");
  panel_file_paths[Porphyromonas_endodontalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Porphyromonas_endodontalis], "data/meta_species/species/Porphyromonas_endodontalis.fa");
  panel_file_paths[Porphyromonas_gingivalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Porphyromonas_gingivalis], "data/meta_species/species/Porphyromonas_gingivalis.fa");
  panel_file_paths[Porphyromonas_uenonis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Porphyromonas_uenonis], "data/meta_species/species/Porphyromonas_uenonis.fa");
  panel_file_paths[Prevotella_amnii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_amnii], "data/meta_species/species/Prevotella_amnii.fa");
  panel_file_paths[Prevotella_bergensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_bergensis], "data/meta_species/species/Prevotella_bergensis.fa");
  panel_file_paths[Prevotella_bivia] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_bivia], "data/meta_species/species/Prevotella_bivia.fa");
  panel_file_paths[Prevotella_bryantii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_bryantii], "data/meta_species/species/Prevotella_bryantii.fa");
  panel_file_paths[Prevotella_buccae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_buccae], "data/meta_species/species/Prevotella_buccae.fa");
  panel_file_paths[Prevotella_buccalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_buccalis], "data/meta_species/species/Prevotella_buccalis.fa");
  panel_file_paths[Prevotella_copri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_copri], "data/meta_species/species/Prevotella_copri.fa");
  panel_file_paths[Prevotella_disiens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_disiens], "data/meta_species/species/Prevotella_disiens.fa");
  panel_file_paths[Prevotella_marshii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_marshii], "data/meta_species/species/Prevotella_marshii.fa");
  panel_file_paths[Prevotella_melaninogenica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_melaninogenica], "data/meta_species/species/Prevotella_melaninogenica.fa");
  panel_file_paths[Prevotella_multiformis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_multiformis], "data/meta_species/species/Prevotella_multiformis.fa");
  panel_file_paths[Prevotella_oralis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_oralis], "data/meta_species/species/Prevotella_oralis.fa");
  panel_file_paths[Prevotella_oris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_oris], "data/meta_species/species/Prevotella_oris.fa");
  panel_file_paths[Prevotella_ruminicola] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_ruminicola], "data/meta_species/species/Prevotella_ruminicola.fa");
  panel_file_paths[Prevotella_salivae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_salivae], "data/meta_species/species/Prevotella_salivae.fa");
  panel_file_paths[Prevotella_tannerae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_tannerae], "data/meta_species/species/Prevotella_tannerae.fa");
  panel_file_paths[Prevotella_timonensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_timonensis], "data/meta_species/species/Prevotella_timonensis.fa");
  panel_file_paths[Prevotella_veroralis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Prevotella_veroralis], "data/meta_species/species/Prevotella_veroralis.fa");
  panel_file_paths[Propionibacterium_acnes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Propionibacterium_acnes], "data/meta_species/species/Propionibacterium_acnes.fa");
  panel_file_paths[Propionibacterium_freudenreichii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Propionibacterium_freudenreichii], "data/meta_species/species/Propionibacterium_freudenreichii.fa");
  panel_file_paths[Proteus_mirabilis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Proteus_mirabilis], "data/meta_species/species/Proteus_mirabilis.fa");
  panel_file_paths[Proteus_penneri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Proteus_penneri], "data/meta_species/species/Proteus_penneri.fa");
  panel_file_paths[Providencia_alcalifaciens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Providencia_alcalifaciens], "data/meta_species/species/Providencia_alcalifaciens.fa");
  panel_file_paths[Providencia_rettgeri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Providencia_rettgeri], "data/meta_species/species/Providencia_rettgeri.fa");
  panel_file_paths[Providencia_rustigianii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Providencia_rustigianii], "data/meta_species/species/Providencia_rustigianii.fa");
  panel_file_paths[Providencia_stuartii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Providencia_stuartii], "data/meta_species/species/Providencia_stuartii.fa");
  panel_file_paths[Pseudoalteromonas_atlantica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudoalteromonas_atlantica], "data/meta_species/species/Pseudoalteromonas_atlantica.fa");
  panel_file_paths[Pseudoalteromonas_haloplanktis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudoalteromonas_haloplanktis], "data/meta_species/species/Pseudoalteromonas_haloplanktis.fa");
  panel_file_paths[Pseudoalteromonas_tunicata] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudoalteromonas_tunicata], "data/meta_species/species/Pseudoalteromonas_tunicata.fa");
  panel_file_paths[Pseudomonas_aeruginosa] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudomonas_aeruginosa], "data/meta_species/species/Pseudomonas_aeruginosa.fa");
  panel_file_paths[Pseudomonas_entomophila] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudomonas_entomophila], "data/meta_species/species/Pseudomonas_entomophila.fa");
  panel_file_paths[Pseudomonas_fluorescens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudomonas_fluorescens], "data/meta_species/species/Pseudomonas_fluorescens.fa");
  panel_file_paths[Pseudomonas_mendocina] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudomonas_mendocina], "data/meta_species/species/Pseudomonas_mendocina.fa");
  panel_file_paths[Pseudomonas_putida] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudomonas_putida], "data/meta_species/species/Pseudomonas_putida.fa");
  panel_file_paths[Pseudomonas_savastanoi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudomonas_savastanoi], "data/meta_species/species/Pseudomonas_savastanoi.fa");
  panel_file_paths[Pseudomonas_stutzeri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudomonas_stutzeri], "data/meta_species/species/Pseudomonas_stutzeri.fa");
  panel_file_paths[Pseudomonas_syringae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pseudomonas_syringae], "data/meta_species/species/Pseudomonas_syringae.fa");
  panel_file_paths[Psychrobacter_arcticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Psychrobacter_arcticus], "data/meta_species/species/Psychrobacter_arcticus.fa");
  panel_file_paths[Psychrobacter_cryohalolentis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Psychrobacter_cryohalolentis], "data/meta_species/species/Psychrobacter_cryohalolentis.fa");
  panel_file_paths[Psychromonas_ingrahamii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Psychromonas_ingrahamii], "data/meta_species/species/Psychromonas_ingrahamii.fa");
  panel_file_paths[Pyrobaculum_aerophilum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pyrobaculum_aerophilum], "data/meta_species/species/Pyrobaculum_aerophilum.fa");
  panel_file_paths[Pyrobaculum_arsenaticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pyrobaculum_arsenaticum], "data/meta_species/species/Pyrobaculum_arsenaticum.fa");
  panel_file_paths[Pyrobaculum_calidifontis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pyrobaculum_calidifontis], "data/meta_species/species/Pyrobaculum_calidifontis.fa");
  panel_file_paths[Pyrobaculum_islandicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pyrobaculum_islandicum], "data/meta_species/species/Pyrobaculum_islandicum.fa");
  panel_file_paths[Pyrococcus_abyssi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pyrococcus_abyssi], "data/meta_species/species/Pyrococcus_abyssi.fa");
  panel_file_paths[Pyrococcus_furiosus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pyrococcus_furiosus], "data/meta_species/species/Pyrococcus_furiosus.fa");
  panel_file_paths[Pyrococcus_horikoshii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Pyrococcus_horikoshii], "data/meta_species/species/Pyrococcus_horikoshii.fa");
  panel_file_paths[Ralstonia_pickettii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ralstonia_pickettii], "data/meta_species/species/Ralstonia_pickettii.fa");
  panel_file_paths[Ralstonia_solanacearum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ralstonia_solanacearum], "data/meta_species/species/Ralstonia_solanacearum.fa");
  panel_file_paths[Rhizobium_etli] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rhizobium_etli], "data/meta_species/species/Rhizobium_etli.fa");
  panel_file_paths[Rhizobium_leguminosarum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rhizobium_leguminosarum], "data/meta_species/species/Rhizobium_leguminosarum.fa");
  panel_file_paths[Rhizobium_radiobacter] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rhizobium_radiobacter], "data/meta_species/species/Rhizobium_radiobacter.fa");
  panel_file_paths[Rhodobacter_capsulatus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rhodobacter_capsulatus], "data/meta_species/species/Rhodobacter_capsulatus.fa");
  panel_file_paths[Rhodobacter_sphaeroides] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rhodobacter_sphaeroides], "data/meta_species/species/Rhodobacter_sphaeroides.fa");
  panel_file_paths[Rhodococcus_equi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rhodococcus_equi], "data/meta_species/species/Rhodococcus_equi.fa");
  panel_file_paths[Rhodococcus_erythropolis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rhodococcus_erythropolis], "data/meta_species/species/Rhodococcus_erythropolis.fa");
  panel_file_paths[Rhodococcus_opacus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rhodococcus_opacus], "data/meta_species/species/Rhodococcus_opacus.fa");
  panel_file_paths[Rickettsia_africae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_africae], "data/meta_species/species/Rickettsia_africae.fa");
  panel_file_paths[Rickettsia_akari] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_akari], "data/meta_species/species/Rickettsia_akari.fa");
  panel_file_paths[Rickettsia_bellii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_bellii], "data/meta_species/species/Rickettsia_bellii.fa");
  panel_file_paths[Rickettsia_canadensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_canadensis], "data/meta_species/species/Rickettsia_canadensis.fa");
  panel_file_paths[Rickettsia_conorii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_conorii], "data/meta_species/species/Rickettsia_conorii.fa");
  panel_file_paths[Rickettsia_endosymbiont] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_endosymbiont], "data/meta_species/species/Rickettsia_endosymbiont.fa");
  panel_file_paths[Rickettsia_felis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_felis], "data/meta_species/species/Rickettsia_felis.fa");
  panel_file_paths[Rickettsia_massiliae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_massiliae], "data/meta_species/species/Rickettsia_massiliae.fa");
  panel_file_paths[Rickettsia_peacockii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_peacockii], "data/meta_species/species/Rickettsia_peacockii.fa");
  panel_file_paths[Rickettsia_prowazekii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_prowazekii], "data/meta_species/species/Rickettsia_prowazekii.fa");
  panel_file_paths[Rickettsia_rickettsii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_rickettsii], "data/meta_species/species/Rickettsia_rickettsii.fa");
  panel_file_paths[Rickettsia_sibirica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_sibirica], "data/meta_species/species/Rickettsia_sibirica.fa");
  panel_file_paths[Rickettsia_typhi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rickettsia_typhi], "data/meta_species/species/Rickettsia_typhi.fa");
  panel_file_paths[Roseburia_intestinalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Roseburia_intestinalis], "data/meta_species/species/Roseburia_intestinalis.fa");
  panel_file_paths[Roseburia_inulinivorans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Roseburia_inulinivorans], "data/meta_species/species/Roseburia_inulinivorans.fa");
  panel_file_paths[Roseiflexus_castenholzii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Roseiflexus_castenholzii], "data/meta_species/species/Roseiflexus_castenholzii.fa");
  panel_file_paths[Roseobacter_denitrificans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Roseobacter_denitrificans], "data/meta_species/species/Roseobacter_denitrificans.fa");
  panel_file_paths[Roseobacter_litoralis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Roseobacter_litoralis], "data/meta_species/species/Roseobacter_litoralis.fa");
  panel_file_paths[Roseovarius_nubinhibens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Roseovarius_nubinhibens], "data/meta_species/species/Roseovarius_nubinhibens.fa");
  panel_file_paths[Rothia_dentocariosa] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rothia_dentocariosa], "data/meta_species/species/Rothia_dentocariosa.fa");
  panel_file_paths[Rothia_mucilaginosa] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Rothia_mucilaginosa], "data/meta_species/species/Rothia_mucilaginosa.fa");
  panel_file_paths[Ruegeria_lacuscaerulensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ruegeria_lacuscaerulensis], "data/meta_species/species/Ruegeria_lacuscaerulensis.fa");
  panel_file_paths[Ruegeria_pomeroyi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ruegeria_pomeroyi], "data/meta_species/species/Ruegeria_pomeroyi.fa");
  panel_file_paths[Ruminococcus_albus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ruminococcus_albus], "data/meta_species/species/Ruminococcus_albus.fa");
  panel_file_paths[Ruminococcus_bromii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ruminococcus_bromii], "data/meta_species/species/Ruminococcus_bromii.fa");
  panel_file_paths[Ruminococcus_flavefaciens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ruminococcus_flavefaciens], "data/meta_species/species/Ruminococcus_flavefaciens.fa");
  panel_file_paths[Ruminococcus_gnavus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ruminococcus_gnavus], "data/meta_species/species/Ruminococcus_gnavus.fa");
  panel_file_paths[Ruminococcus_lactaris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ruminococcus_lactaris], "data/meta_species/species/Ruminococcus_lactaris.fa");
  panel_file_paths[Ruminococcus_obeum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ruminococcus_obeum], "data/meta_species/species/Ruminococcus_obeum.fa");
  panel_file_paths[Ruminococcus_torques] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ruminococcus_torques], "data/meta_species/species/Ruminococcus_torques.fa");
  panel_file_paths[Salinispora_arenicola] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Salinispora_arenicola], "data/meta_species/species/Salinispora_arenicola.fa");
  panel_file_paths[Salinispora_tropica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Salinispora_tropica], "data/meta_species/species/Salinispora_tropica.fa");
  panel_file_paths[Salmonella_enterica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Salmonella_enterica], "data/meta_species/species/Salmonella_enterica.fa");
  panel_file_paths[Salmonella_typhimurium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Salmonella_typhimurium], "data/meta_species/species/Salmonella_typhimurium.fa");
  panel_file_paths[Segniliparus_rotundus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Segniliparus_rotundus], "data/meta_species/species/Segniliparus_rotundus.fa");
  panel_file_paths[Segniliparus_rugosus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Segniliparus_rugosus], "data/meta_species/species/Segniliparus_rugosus.fa");
  panel_file_paths[Selenomonas_artemidis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Selenomonas_artemidis], "data/meta_species/species/Selenomonas_artemidis.fa");
  panel_file_paths[Selenomonas_flueggei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Selenomonas_flueggei], "data/meta_species/species/Selenomonas_flueggei.fa");
  panel_file_paths[Selenomonas_noxia] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Selenomonas_noxia], "data/meta_species/species/Selenomonas_noxia.fa");
  panel_file_paths[Selenomonas_sputigena] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Selenomonas_sputigena], "data/meta_species/species/Selenomonas_sputigena.fa");
  panel_file_paths[Serratia_odorifera] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Serratia_odorifera], "data/meta_species/species/Serratia_odorifera.fa");
  panel_file_paths[Serratia_proteamaculans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Serratia_proteamaculans], "data/meta_species/species/Serratia_proteamaculans.fa");
  panel_file_paths[Shewanella_amazonensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_amazonensis], "data/meta_species/species/Shewanella_amazonensis.fa");
  panel_file_paths[Shewanella_baltica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_baltica], "data/meta_species/species/Shewanella_baltica.fa");
  panel_file_paths[Shewanella_benthica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_benthica], "data/meta_species/species/Shewanella_benthica.fa");
  panel_file_paths[Shewanella_denitrificans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_denitrificans], "data/meta_species/species/Shewanella_denitrificans.fa");
  panel_file_paths[Shewanella_frigidimarina] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_frigidimarina], "data/meta_species/species/Shewanella_frigidimarina.fa");
  panel_file_paths[Shewanella_halifaxensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_halifaxensis], "data/meta_species/species/Shewanella_halifaxensis.fa");
  panel_file_paths[Shewanella_loihica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_loihica], "data/meta_species/species/Shewanella_loihica.fa");
  panel_file_paths[Shewanella_oneidensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_oneidensis], "data/meta_species/species/Shewanella_oneidensis.fa");
  panel_file_paths[Shewanella_pealeana] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_pealeana], "data/meta_species/species/Shewanella_pealeana.fa");
  panel_file_paths[Shewanella_piezotolerans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_piezotolerans], "data/meta_species/species/Shewanella_piezotolerans.fa");
  panel_file_paths[Shewanella_putrefaciens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_putrefaciens], "data/meta_species/species/Shewanella_putrefaciens.fa");
  panel_file_paths[Shewanella_sediminis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_sediminis], "data/meta_species/species/Shewanella_sediminis.fa");
  panel_file_paths[Shewanella_violacea] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_violacea], "data/meta_species/species/Shewanella_violacea.fa");
  panel_file_paths[Shewanella_woodyi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shewanella_woodyi], "data/meta_species/species/Shewanella_woodyi.fa");
  panel_file_paths[Shigella_boydii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shigella_boydii], "data/meta_species/species/Shigella_boydii.fa");
  panel_file_paths[Shigella_dysenteriae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shigella_dysenteriae], "data/meta_species/species/Shigella_dysenteriae.fa");
  panel_file_paths[Shigella_flexneri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shigella_flexneri], "data/meta_species/species/Shigella_flexneri.fa");
  panel_file_paths[Shigella_sonnei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shigella_sonnei], "data/meta_species/species/Shigella_sonnei.fa");
  panel_file_paths[Slackia_exigua] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Slackia_exigua], "data/meta_species/species/Slackia_exigua.fa");
  panel_file_paths[Slackia_heliotrinireducens] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Slackia_heliotrinireducens], "data/meta_species/species/Slackia_heliotrinireducens.fa");
  panel_file_paths[Sphingobium_chlorophenolicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sphingobium_chlorophenolicum], "data/meta_species/species/Sphingobium_chlorophenolicum.fa");
  panel_file_paths[Sphingobium_japonicum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sphingobium_japonicum], "data/meta_species/species/Sphingobium_japonicum.fa");
  panel_file_paths[Sphingomonas_wittichii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sphingomonas_wittichii], "data/meta_species/species/Sphingomonas_wittichii.fa");
  panel_file_paths[Spirochaeta_smaragdinae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Spirochaeta_smaragdinae], "data/meta_species/species/Spirochaeta_smaragdinae.fa");
  panel_file_paths[Spirochaeta_thermophila] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Spirochaeta_thermophila], "data/meta_species/species/Spirochaeta_thermophila.fa");
  panel_file_paths[Saureus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Saureus], "data/meta_species/species/Staphylococcus_aureus.fa");
  panel_file_paths[Scapitis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Scapitis], "data/meta_species/species/Staphylococcus_capitis.fa");
  panel_file_paths[Scarnosus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Scarnosus], "data/meta_species/species/Staphylococcus_carnosus.fa");
  panel_file_paths[Sepidermidis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sepidermidis], "data/meta_species/species/Staphylococcus_epidermidis.fa");
  panel_file_paths[Shaemolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shaemolyticus], "data/meta_species/species/Staphylococcus_haemolyticus.fa");
  panel_file_paths[Shominis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Shominis], "data/meta_species/species/Staphylococcus_hominis.fa");
  panel_file_paths[Slugdunensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Slugdunensis], "data/meta_species/species/Staphylococcus_lugdunensis.fa");
  panel_file_paths[Spseudintermedius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Spseudintermedius], "data/meta_species/species/Staphylococcus_pseudintermedius.fa");
  panel_file_paths[Ssaprophyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ssaprophyticus], "data/meta_species/species/Staphylococcus_saprophyticus.fa");
  panel_file_paths[Swarneri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Swarneri], "data/meta_species/species/Staphylococcus_warneri.fa");
  panel_file_paths[Staphylothermus_hellenicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Staphylothermus_hellenicus], "data/meta_species/species/Staphylothermus_hellenicus.fa");
  panel_file_paths[Staphylothermus_marinus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Staphylothermus_marinus], "data/meta_species/species/Staphylothermus_marinus.fa");
  panel_file_paths[Stenotrophomonas_maltophilia] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Stenotrophomonas_maltophilia], "data/meta_species/species/Stenotrophomonas_maltophilia.fa");
  panel_file_paths[Streptococcus_agalactiae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_agalactiae], "data/meta_species/species/Streptococcus_agalactiae.fa");
  panel_file_paths[Streptococcus_anginosus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_anginosus], "data/meta_species/species/Streptococcus_anginosus.fa");
  panel_file_paths[Streptococcus_australis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_australis], "data/meta_species/species/Streptococcus_australis.fa");
  panel_file_paths[Streptococcus_bovis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_bovis], "data/meta_species/species/Streptococcus_bovis.fa");
  panel_file_paths[Streptococcus_cristatus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_cristatus], "data/meta_species/species/Streptococcus_cristatus.fa");
  panel_file_paths[Streptococcus_downei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_downei], "data/meta_species/species/Streptococcus_downei.fa");
  panel_file_paths[Streptococcus_dysgalactiae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_dysgalactiae], "data/meta_species/species/Streptococcus_dysgalactiae.fa");
  panel_file_paths[Streptococcus_equi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_equi], "data/meta_species/species/Streptococcus_equi.fa");
  panel_file_paths[Streptococcus_equinus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_equinus], "data/meta_species/species/Streptococcus_equinus.fa");
  panel_file_paths[Streptococcus_gallolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_gallolyticus], "data/meta_species/species/Streptococcus_gallolyticus.fa");
  panel_file_paths[Streptococcus_gordonii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_gordonii], "data/meta_species/species/Streptococcus_gordonii.fa");
  panel_file_paths[Streptococcus_infantarius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_infantarius], "data/meta_species/species/Streptococcus_infantarius.fa");
  panel_file_paths[Streptococcus_infantis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_infantis], "data/meta_species/species/Streptococcus_infantis.fa");
  panel_file_paths[Streptococcus_mitis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_mitis], "data/meta_species/species/Streptococcus_mitis.fa");
  panel_file_paths[Streptococcus_mutans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_mutans], "data/meta_species/species/Streptococcus_mutans.fa");
  panel_file_paths[Streptococcus_oralis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_oralis], "data/meta_species/species/Streptococcus_oralis.fa");
  panel_file_paths[Streptococcus_parasanguinis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_parasanguinis], "data/meta_species/species/Streptococcus_parasanguinis.fa");
  panel_file_paths[Streptococcus_peroris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_peroris], "data/meta_species/species/Streptococcus_peroris.fa");
  panel_file_paths[Streptococcus_pneumoniae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_pneumoniae], "data/meta_species/species/Streptococcus_pneumoniae.fa");
  panel_file_paths[Streptococcus_pseudoporcinus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_pseudoporcinus], "data/meta_species/species/Streptococcus_pseudoporcinus.fa");
  panel_file_paths[Streptococcus_pyogenes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_pyogenes], "data/meta_species/species/Streptococcus_pyogenes.fa");
  panel_file_paths[Streptococcus_salivarius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_salivarius], "data/meta_species/species/Streptococcus_salivarius.fa");
  panel_file_paths[Streptococcus_sanguinis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_sanguinis], "data/meta_species/species/Streptococcus_sanguinis.fa");
  panel_file_paths[Streptococcus_suis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_suis], "data/meta_species/species/Streptococcus_suis.fa");
  panel_file_paths[Streptococcus_thermophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_thermophilus], "data/meta_species/species/Streptococcus_thermophilus.fa");
  panel_file_paths[Streptococcus_uberis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_uberis], "data/meta_species/species/Streptococcus_uberis.fa");
  panel_file_paths[Streptococcus_vestibularis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptococcus_vestibularis], "data/meta_species/species/Streptococcus_vestibularis.fa");
  panel_file_paths[Streptomyces_albus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_albus], "data/meta_species/species/Streptomyces_albus.fa");
  panel_file_paths[Streptomyces_avermitilis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_avermitilis], "data/meta_species/species/Streptomyces_avermitilis.fa");
  panel_file_paths[Streptomyces_bingchenggensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_bingchenggensis], "data/meta_species/species/Streptomyces_bingchenggensis.fa");
  panel_file_paths[Streptomyces_clavuligerus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_clavuligerus], "data/meta_species/species/Streptomyces_clavuligerus.fa");
  panel_file_paths[Streptomyces_coelicolor] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_coelicolor], "data/meta_species/species/Streptomyces_coelicolor.fa");
  panel_file_paths[Streptomyces_ghanaensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_ghanaensis], "data/meta_species/species/Streptomyces_ghanaensis.fa");
  panel_file_paths[Streptomyces_griseoflavus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_griseoflavus], "data/meta_species/species/Streptomyces_griseoflavus.fa");
  panel_file_paths[Streptomyces_griseus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_griseus], "data/meta_species/species/Streptomyces_griseus.fa");
  panel_file_paths[Streptomyces_hygroscopicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_hygroscopicus], "data/meta_species/species/Streptomyces_hygroscopicus.fa");
  panel_file_paths[Streptomyces_lividans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_lividans], "data/meta_species/species/Streptomyces_lividans.fa");
  panel_file_paths[Streptomyces_pristinaespiralis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_pristinaespiralis], "data/meta_species/species/Streptomyces_pristinaespiralis.fa");
  panel_file_paths[Streptomyces_roseosporus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_roseosporus], "data/meta_species/species/Streptomyces_roseosporus.fa");
  panel_file_paths[Streptomyces_scabiei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_scabiei], "data/meta_species/species/Streptomyces_scabiei.fa");
  panel_file_paths[Streptomyces_sviceus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_sviceus], "data/meta_species/species/Streptomyces_sviceus.fa");
  panel_file_paths[Streptomyces_violaceusniger] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_violaceusniger], "data/meta_species/species/Streptomyces_violaceusniger.fa");
  panel_file_paths[Streptomyces_viridochromogenes] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Streptomyces_viridochromogenes], "data/meta_species/species/Streptomyces_viridochromogenes.fa");
  panel_file_paths[Sulfolobus_acidocaldarius] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sulfolobus_acidocaldarius], "data/meta_species/species/Sulfolobus_acidocaldarius.fa");
  panel_file_paths[Sulfolobus_islandicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sulfolobus_islandicus], "data/meta_species/species/Sulfolobus_islandicus.fa");
  panel_file_paths[Sulfolobus_solfataricus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sulfolobus_solfataricus], "data/meta_species/species/Sulfolobus_solfataricus.fa");
  panel_file_paths[Sulfolobus_tokodaii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sulfolobus_tokodaii], "data/meta_species/species/Sulfolobus_tokodaii.fa");
  panel_file_paths[Sulfurihydrogenibium_azorense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sulfurihydrogenibium_azorense], "data/meta_species/species/Sulfurihydrogenibium_azorense.fa");
  panel_file_paths[Sulfurihydrogenibium_yellowstonense] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sulfurihydrogenibium_yellowstonense], "data/meta_species/species/Sulfurihydrogenibium_yellowstonense.fa");
  panel_file_paths[Sulfurimonas_autotrophica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sulfurimonas_autotrophica], "data/meta_species/species/Sulfurimonas_autotrophica.fa");
  panel_file_paths[Sulfurimonas_denitrificans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Sulfurimonas_denitrificans], "data/meta_species/species/Sulfurimonas_denitrificans.fa");
  panel_file_paths[Synechococcus_elongatus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Synechococcus_elongatus], "data/meta_species/species/Synechococcus_elongatus.fa");
  panel_file_paths[Thermaerobacter_marianensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermaerobacter_marianensis], "data/meta_species/species/Thermaerobacter_marianensis.fa");
  panel_file_paths[Thermaerobacter_subterraneus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermaerobacter_subterraneus], "data/meta_species/species/Thermaerobacter_subterraneus.fa");
  panel_file_paths[Thermoanaerobacter_brockii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoanaerobacter_brockii], "data/meta_species/species/Thermoanaerobacter_brockii.fa");
  panel_file_paths[Thermoanaerobacter_ethanolicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoanaerobacter_ethanolicus], "data/meta_species/species/Thermoanaerobacter_ethanolicus.fa");
  panel_file_paths[Thermoanaerobacter_italicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoanaerobacter_italicus], "data/meta_species/species/Thermoanaerobacter_italicus.fa");
  panel_file_paths[Thermoanaerobacterium_thermosaccharolyticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoanaerobacterium_thermosaccharolyticum], "data/meta_species/species/Thermoanaerobacterium_thermosaccharolyticum.fa");
  panel_file_paths[Thermoanaerobacterium_xylanolyticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoanaerobacterium_xylanolyticum], "data/meta_species/species/Thermoanaerobacterium_xylanolyticum.fa");
  panel_file_paths[Thermoanaerobacter_mathranii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoanaerobacter_mathranii], "data/meta_species/species/Thermoanaerobacter_mathranii.fa");
  panel_file_paths[Thermoanaerobacter_pseudethanolicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoanaerobacter_pseudethanolicus], "data/meta_species/species/Thermoanaerobacter_pseudethanolicus.fa");
  panel_file_paths[Thermoanaerobacter_wiegelii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoanaerobacter_wiegelii], "data/meta_species/species/Thermoanaerobacter_wiegelii.fa");
  panel_file_paths[Thermococcus_barophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermococcus_barophilus], "data/meta_species/species/Thermococcus_barophilus.fa");
  panel_file_paths[Thermococcus_gammatolerans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermococcus_gammatolerans], "data/meta_species/species/Thermococcus_gammatolerans.fa");
  panel_file_paths[Thermococcus_kodakaraensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermococcus_kodakaraensis], "data/meta_species/species/Thermococcus_kodakaraensis.fa");
  panel_file_paths[Thermococcus_onnurineus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermococcus_onnurineus], "data/meta_species/species/Thermococcus_onnurineus.fa");
  panel_file_paths[Thermococcus_sibiricus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermococcus_sibiricus], "data/meta_species/species/Thermococcus_sibiricus.fa");
  panel_file_paths[Thermoplasma_acidophilum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoplasma_acidophilum], "data/meta_species/species/Thermoplasma_acidophilum.fa");
  panel_file_paths[Thermoplasma_volcanium] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermoplasma_volcanium], "data/meta_species/species/Thermoplasma_volcanium.fa");
  panel_file_paths[Thermosipho_africanus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermosipho_africanus], "data/meta_species/species/Thermosipho_africanus.fa");
  panel_file_paths[Thermosipho_melanesiensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermosipho_melanesiensis], "data/meta_species/species/Thermosipho_melanesiensis.fa");
  panel_file_paths[Thermotoga_lettingae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermotoga_lettingae], "data/meta_species/species/Thermotoga_lettingae.fa");
  panel_file_paths[Thermotoga_maritima] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermotoga_maritima], "data/meta_species/species/Thermotoga_maritima.fa");
  panel_file_paths[Thermotoga_naphthophila] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermotoga_naphthophila], "data/meta_species/species/Thermotoga_naphthophila.fa");
  panel_file_paths[Thermotoga_neapolitana] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermotoga_neapolitana], "data/meta_species/species/Thermotoga_neapolitana.fa");
  panel_file_paths[Thermotoga_petrophila] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermotoga_petrophila], "data/meta_species/species/Thermotoga_petrophila.fa");
  panel_file_paths[Thermus_aquaticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermus_aquaticus], "data/meta_species/species/Thermus_aquaticus.fa");
  panel_file_paths[Thermus_scotoductus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermus_scotoductus], "data/meta_species/species/Thermus_scotoductus.fa");
  panel_file_paths[Thermus_thermophilus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Thermus_thermophilus], "data/meta_species/species/Thermus_thermophilus.fa");
  panel_file_paths[Treponema_denticola] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Treponema_denticola], "data/meta_species/species/Treponema_denticola.fa");
  panel_file_paths[Treponema_pallidum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Treponema_pallidum], "data/meta_species/species/Treponema_pallidum.fa");
  panel_file_paths[Treponema_phagedenis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Treponema_phagedenis], "data/meta_species/species/Treponema_phagedenis.fa");
  panel_file_paths[Treponema_vincentii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Treponema_vincentii], "data/meta_species/species/Treponema_vincentii.fa");
  panel_file_paths[Ureaplasma_parvum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ureaplasma_parvum], "data/meta_species/species/Ureaplasma_parvum.fa");
  panel_file_paths[Ureaplasma_urealyticum] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Ureaplasma_urealyticum], "data/meta_species/species/Ureaplasma_urealyticum.fa");
  panel_file_paths[Veillonella_atypica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Veillonella_atypica], "data/meta_species/species/Veillonella_atypica.fa");
  panel_file_paths[Veillonella_dispar] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Veillonella_dispar], "data/meta_species/species/Veillonella_dispar.fa");
  panel_file_paths[Veillonella_parvula] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Veillonella_parvula], "data/meta_species/species/Veillonella_parvula.fa");
  panel_file_paths[Vibrio_alginolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_alginolyticus], "data/meta_species/species/Vibrio_alginolyticus.fa");
  panel_file_paths[Vibrio_brasiliensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_brasiliensis], "data/meta_species/species/Vibrio_brasiliensis.fa");
  panel_file_paths[Vibrio_campbellii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_campbellii], "data/meta_species/species/Vibrio_campbellii.fa");
  panel_file_paths[Vibrio_caribbenthicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_caribbenthicus], "data/meta_species/species/Vibrio_caribbenthicus.fa");
  panel_file_paths[Vibrio_cholerae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_cholerae], "data/meta_species/species/Vibrio_cholerae.fa");
  panel_file_paths[Vibrio_coralliilyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_coralliilyticus], "data/meta_species/species/Vibrio_coralliilyticus.fa");
  panel_file_paths[Vibrio_furnissii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_furnissii], "data/meta_species/species/Vibrio_furnissii.fa");
  panel_file_paths[Vibrio_harveyi] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_harveyi], "data/meta_species/species/Vibrio_harveyi.fa");
  panel_file_paths[Vibrio_metschnikovii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_metschnikovii], "data/meta_species/species/Vibrio_metschnikovii.fa");
  panel_file_paths[Vibrio_mimicus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_mimicus], "data/meta_species/species/Vibrio_mimicus.fa");
  panel_file_paths[Vibrio_orientalis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_orientalis], "data/meta_species/species/Vibrio_orientalis.fa");
  panel_file_paths[Vibrio_parahaemolyticus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_parahaemolyticus], "data/meta_species/species/Vibrio_parahaemolyticus.fa");
  panel_file_paths[Vibrio_shilonii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_shilonii], "data/meta_species/species/Vibrio_shilonii.fa");
  panel_file_paths[Vibrio_sinaloensis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_sinaloensis], "data/meta_species/species/Vibrio_sinaloensis.fa");
  panel_file_paths[Vibrio_splendidus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_splendidus], "data/meta_species/species/Vibrio_splendidus.fa");
  panel_file_paths[Vibrio_vulnificus] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vibrio_vulnificus], "data/meta_species/species/Vibrio_vulnificus.fa");
  panel_file_paths[Vulcanisaeta_distributa] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vulcanisaeta_distributa], "data/meta_species/species/Vulcanisaeta_distributa.fa");
  panel_file_paths[Vulcanisaeta_moutnovskia] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Vulcanisaeta_moutnovskia], "data/meta_species/species/Vulcanisaeta_moutnovskia.fa");
  panel_file_paths[Wolbachia_endosymbiont] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Wolbachia_endosymbiont], "data/meta_species/species/Wolbachia_endosymbiont.fa");
  panel_file_paths[Wolbachia_pipientis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Wolbachia_pipientis], "data/meta_species/species/Wolbachia_pipientis.fa");
  panel_file_paths[Xanthomonas_albilineans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Xanthomonas_albilineans], "data/meta_species/species/Xanthomonas_albilineans.fa");
  panel_file_paths[Xanthomonas_axonopodis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Xanthomonas_axonopodis], "data/meta_species/species/Xanthomonas_axonopodis.fa");
  panel_file_paths[Xanthomonas_campestris] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Xanthomonas_campestris], "data/meta_species/species/Xanthomonas_campestris.fa");
  panel_file_paths[Xanthomonas_fuscans] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Xanthomonas_fuscans], "data/meta_species/species/Xanthomonas_fuscans.fa");
  panel_file_paths[Xanthomonas_oryzae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Xanthomonas_oryzae], "data/meta_species/species/Xanthomonas_oryzae.fa");
  panel_file_paths[Xenorhabdus_bovienii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Xenorhabdus_bovienii], "data/meta_species/species/Xenorhabdus_bovienii.fa");
  panel_file_paths[Xenorhabdus_nematophila] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Xenorhabdus_nematophila], "data/meta_species/species/Xenorhabdus_nematophila.fa");
  panel_file_paths[Yersinia_aldovae] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_aldovae], "data/meta_species/species/Yersinia_aldovae.fa");
  panel_file_paths[Yersinia_bercovieri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_bercovieri], "data/meta_species/species/Yersinia_bercovieri.fa");
  panel_file_paths[Yersinia_enterocolitica] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_enterocolitica], "data/meta_species/species/Yersinia_enterocolitica.fa");
  panel_file_paths[Yersinia_frederiksenii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_frederiksenii], "data/meta_species/species/Yersinia_frederiksenii.fa");
  panel_file_paths[Yersinia_intermedia] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_intermedia], "data/meta_species/species/Yersinia_intermedia.fa");
  panel_file_paths[Yersinia_kristensenii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_kristensenii], "data/meta_species/species/Yersinia_kristensenii.fa");
  panel_file_paths[Yersinia_mollaretii] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_mollaretii], "data/meta_species/species/Yersinia_mollaretii.fa");
  panel_file_paths[Yersinia_pestis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_pestis], "data/meta_species/species/Yersinia_pestis.fa");
  panel_file_paths[Yersinia_pseudotuberculosis] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_pseudotuberculosis], "data/meta_species/species/Yersinia_pseudotuberculosis.fa");
  panel_file_paths[Yersinia_rohdei] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_rohdei], "data/meta_species/species/Yersinia_rohdei.fa");
  panel_file_paths[Yersinia_ruckeri] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[Yersinia_ruckeri], "data/meta_species/species/Yersinia_ruckeri.fa");
}

void phylo_group_threshold(int* thresholds){
  thresholds[CoagPos] = 90;
  thresholds[CoagNeg] = 25;

}

void cat_threshold(int* thresholds){
  thresholds[0] = 20;
}
void load_all_species_thresholds(int* thresholds){
  int j;
  for(j = 0; j < NUM_SPECIES; j++) {
    thresholds[j] = 90;
  }   
  thresholds[Saureus] = 90;
  thresholds[Sepidermidis] = 30;
  thresholds[Shaemolyticus] = 30;
  thresholds[Sother] = 10;
}

boolean* create_staph_mask()
{
  boolean* mask= create_mask(false);
  mask[Saureus] = true;
  mask[Scapitis] = true;
  mask[Scarnosus] = true;
  mask[Sepidermidis] = true;
  mask[Shaemolyticus] = true;
  mask[Shominis] = true;
  mask[Slugdunensis] = true;
  mask[Spseudintermedius] = true;
  mask[Ssaprophyticus] = true;
  mask[Swarneri] = true;
  return (mask);
}

boolean* create_non_aureus_mask()
{
  boolean* mask= create_staph_mask();
  mask[Saureus] = false;
  return (mask);
}

boolean non_aureus_panels_are_present(SpeciesInfo* species_info){
  boolean* mask = create_non_aureus_mask();
  boolean non_aureus_species_panels_are_present = panels_are_present(species_info->species_covg_info,mask);
  return (non_aureus_species_panels_are_present);
}

boolean no_non_aureus_panels_are_present(SpeciesInfo* species_info){
  return (!non_aureus_panels_are_present(species_info));
}

boolean staphylococcus_is_present(SpeciesInfo* species_info){
  boolean* mask = create_staph_mask();
  boolean staph_species_panels_are_present = panels_are_present(species_info->species_covg_info,mask);
  return (staph_species_panels_are_present);
}

Species get_best_staph_species(SpeciesInfo* species_info ){
  boolean* mask = create_staph_mask();
  int species_enum = get_best_hit(species_info->species_covg_info,mask);
  Species species = species_enum;
  return (species);
}

Species get_best_species(SpeciesInfo* species_info ){
  boolean* mask = create_mask(true);
  int species_enum = get_best_hit(species_info->species_covg_info,mask);
  Species species = species_enum;
  return (species);
}

Species get_best_non_aureus_species(SpeciesInfo* species_info ){
  boolean* mask = create_non_aureus_mask();
  int species_enum  = get_best_hit(species_info->species_covg_info,mask);
  Species species = species_enum;
  return (species);
}


boolean is_aureus_present(SpeciesInfo* species_info)
{
  return (species_info->species_covg_info->present[Saureus]);
}

boolean is_non_aureus_staph_present(SpeciesInfo* species_info)
{
  boolean is_epi_present = species_info->species_covg_info->present[Sepidermidis];
  boolean is_haem_present = species_info->species_covg_info->present[Shaemolyticus];
  boolean is_sother_present = species_info->species_covg_info->present[Sother];
  if (is_epi_present || is_haem_present  || is_sother_present){
    return (true);
  }
  else{
    return (false);
  }
}

void update_phylo_group_presence_and_coverage_from_species(SpeciesInfo* species_info){
  if (non_aureus_panels_are_present(species_info)){
    if (! species_info->phylo_group_covg_info->present[CoagNeg]){
      species_info->phylo_group_covg_info->present[CoagNeg] = true;
      species_info->phylo_group_covg_info->num_panels_present = species_info->phylo_group_covg_info->num_panels_present + 1;
    }
    Species best_staph_species = get_best_non_aureus_species(species_info);
    species_info->phylo_group_covg_info->percentage_coverage[CoagNeg] = max(species_info->phylo_group_covg_info->percentage_coverage[CoagNeg] , species_info->species_covg_info->percentage_coverage[best_staph_species] );
    species_info->phylo_group_covg_info->median_coverage[CoagNeg] = max(species_info->phylo_group_covg_info->median_coverage[CoagNeg] , species_info->species_covg_info->median_coverage[best_staph_species] );
  }
}

SpeciesInfo* get_species_info(dBGraph *db_graph,int max_branch_len, 
                            StrBuf* install_dir,int expected_covg,
                            int ignore_first,int ignore_last)

{
  StrBuf* phylo_group_file_paths[NUM_COMPLEX];
  load_all_phylo_group_file_paths(phylo_group_file_paths,install_dir);
  StrBuf* species_file_paths[NUM_SPECIES];
  load_all_species_file_paths(species_file_paths,install_dir);


  CovgInfo* phylo_group_covg_info = get_coverage_info(db_graph,
                                                  phylo_group_file_paths,
                                                  max_branch_len,NUM_COMPLEX,
                                                  ignore_first,ignore_last,
                                                  phylo_group_threshold);
  CovgInfo* species_covg_info = get_coverage_info(db_graph,
                                                  species_file_paths,
                                                  max_branch_len,NUM_SPECIES,
                                                  ignore_first,ignore_last,
                                                  load_all_species_thresholds);

  StrBuf* cat_file_paths[1];
  load_all_cat_file_paths(cat_file_paths,install_dir);    
  CovgInfo* cat_covg_info = get_coverage_info(db_graph,
                                            cat_file_paths,
                                            max_branch_len,1,
                                            ignore_first,ignore_last,
                                            cat_threshold);   



  SpeciesInfo* species_info=(SpeciesInfo *)malloc(sizeof(SpeciesInfo)); 
  species_info->phylo_group_covg_info = phylo_group_covg_info;
  species_info->species_covg_info = species_covg_info;
  species_info->other_covg_info = cat_covg_info;

  update_phylo_group_presence_and_coverage_from_species(species_info);

  return species_info;
}



void print_json_aureus(SpeciesInfo* species_info, boolean last){
    print_json_called_variant_item( get_char_name_of_species_enum (Saureus) ,species_info->species_covg_info->median_coverage[Saureus], last);
}

void print_json_best_hit_non_aureus(SpeciesInfo* species_info){
  if (no_non_aureus_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Species", -1 , true);
  }
  else{
  Species staph_species = get_best_non_aureus_species(species_info);
  print_json_called_variant_item( get_char_name_of_species_enum(staph_species), species_info->species_covg_info->median_coverage[staph_species], true);    
  }
}

void print_json_best_hit(SpeciesInfo* species_info){
  Species species = get_best_species(species_info);
  print_json_called_variant_item( get_char_name_of_species_enum(species), species_info->species_covg_info->median_coverage[species], true);    
}  

void print_json_aureus_and_best_hit_non_aureus(SpeciesInfo* species_info){
  if (is_aureus_present(species_info)){
    print_json_aureus(species_info,false);
  }
  else
  {  
    print_json_called_variant_item( "Unknown Species", -1 , false);
  }
  if (no_non_aureus_panels_are_present(species_info)){
    print_json_called_variant_item( "Unknown Species", -1 , true);
  }
  else
  {  
    Species staph_species = get_best_non_aureus_species(species_info);  
    print_json_called_variant_item( get_char_name_of_species_enum(staph_species), species_info->species_covg_info->median_coverage[staph_species], true);    
  }
}



boolean catalayse_exists_in_sample(SpeciesInfo* species_info)

{    
    if (species_info->other_covg_info->percentage_coverage[0] > 20)
    {
      return true;
    }else
    {
      return false;
    }
}

int get_coverage_on_catalayse(SpeciesInfo* species_info)
{    
  return(species_info->other_covg_info->median_coverage[0]);
}

void load_all_cat_file_paths(StrBuf** panel_file_paths , StrBuf* install_dir )
{
  panel_file_paths[0] = strbuf_create(install_dir->buff);
  strbuf_append_str(panel_file_paths[0], "data/staph/species/coag_neg.fasta" ); 
}

void print_json_phylo_group(SpeciesInfo* species_info){
    CovgInfo* covg_info =species_info->phylo_group_covg_info;    
    int num_panels_present = covg_info->num_panels_present;
    print_json_phylo_group_start();
    if (num_panels_present > 0){
      print_json_indiv_phylo(covg_info,get_ith_phylo_group_name);
    }
    else
    {
      if (catalayse_exists_in_sample(species_info)){
        print_json_called_variant_item( "Coagulase-Negative Staphylococcus", get_coverage_on_catalayse(species_info), true);
      }
      else{
        print_json_called_variant_item( "Non Staphylococcus", -1, true);
      }
    }
    print_json_phylo_group_end();  
}


void print_json_species(SpeciesInfo* species_info){
    Species aureus_is_present = is_aureus_present(species_info);
    Species non_aureus_staph_is_present = is_non_aureus_staph_present(species_info);
    print_json_species_start();
    if (aureus_is_present && non_aureus_staph_is_present){
      print_json_aureus_and_best_hit_non_aureus(species_info);
    }
    else if (aureus_is_present){
      print_json_aureus(species_info,true);
    }
    else if (non_aureus_staph_is_present){
      print_json_best_hit_non_aureus(species_info);
    }
    else
    {
      int num_panels_present = species_info->species_covg_info->num_panels_present;
      if (num_panels_present > 0 ){
        print_json_best_hit( species_info);
      }else{
        print_json_called_variant_item( "Unknown Species", -1, true);
      }
      
    }    
    print_json_species_end();  
}

void print_json_lineage(SpeciesInfo* species_info){
    print_json_lineage_start();
    print_json_called_variant_item( "N/A", -1, true);
    print_json_lineage_end(); 
}
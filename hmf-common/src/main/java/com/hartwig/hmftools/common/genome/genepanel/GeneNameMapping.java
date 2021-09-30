package com.hartwig.hmftools.common.genome.genepanel;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

// a resource to map gene names between v37 and v38
// if a gene is not present, then it's name has not changed
public class GeneNameMapping
{
    public static final Set<String> IG_GENES = Sets.newHashSet("IGH", "IGK", "IGL");

    private final Map<String, String> m37to38Map;
    private final Map<String, String> m38to37Map;
    private final Set<String> mUnchangedGenes;

    private static final String DELIM = "\t";

    @NotNull
    public static GeneNameMapping loadFromEmbeddedResource()
    {
        final Map<String, String> v37Map = Maps.newHashMap();
        final Map<String, String> v38Map = Maps.newHashMap();
        final Set<String> unchangedGenes = Sets.newHashSet();

        final InputStream inputStream = GeneNameMapping.class.getResourceAsStream("/ensembl/gene_name_mapping_37_38.tsv");
        List<String> mappingLines = new BufferedReader(new InputStreamReader(inputStream)).lines().collect(Collectors.toList());
        mappingLines.remove(0);

        for(String mapping : mappingLines)
        {
            String[] values = mapping.split(DELIM);
            String geneName37 = values[0];
            String geneName38 = values[1];

            if(geneName37.equals(geneName38))
            {
                unchangedGenes.add(geneName37);
            }
            else
            {
                v37Map.put(geneName37, geneName38);
                v38Map.put(geneName38, geneName37);
            }
        }

        return new GeneNameMapping(v37Map, v38Map, unchangedGenes);
    }

    public GeneNameMapping(final Map<String,String> v37Map, final Map<String,String> v38Map, final Set<String> unchangedGenes)
    {
        m37to38Map = v37Map;
        m38to37Map = v38Map;
        mUnchangedGenes = unchangedGenes;
    }

    public boolean isValidV38Gene(@NotNull final String v38GeneId)
    {
        return mUnchangedGenes.contains(v38GeneId)
            || (m38to37Map.containsKey(v38GeneId) && !m38to37Map.get(v38GeneId).equals("NA"))
            || IG_GENES.contains(v38GeneId);
    }

    public boolean isValidV37Gene(@NotNull final String v37GeneId)
    {
        return mUnchangedGenes.contains(v37GeneId) || m37to38Map.containsKey(v37GeneId) || IG_GENES.contains(v37GeneId);
    }

    @NotNull
    public String v37Gene(@NotNull final String geneNamev38)
    {
        if(IG_GENES.contains(geneNamev38) || mUnchangedGenes.contains(geneNamev38))
            return geneNamev38;

        String result = m38to37Map.get(geneNamev38);
        return result != null ? result : geneNamev38;
    }

    @NotNull
    public String v38Gene(@NotNull final String geneNamev37)
    {
        if(IG_GENES.contains(geneNamev37) || mUnchangedGenes.contains(geneNamev37))
            return geneNamev37;

        String result = m37to38Map.get(geneNamev37);
        return result != null ? result : geneNamev37;
    }
}

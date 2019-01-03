package com.hartwig.hmftools.common.variant.structural.annotation;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

public class EnsemblGeneData
{
    public final String GeneId;
    public final String GeneName;
    public final String Chromosome;
    public final byte Strand;
    public final long GeneStart;
    public final long GeneEnd;
    public final List<Integer> EntrezIds;
    public final String KaryotypeBand;
    public final List<String> Synonyms;

    public EnsemblGeneData(String geneId, String geneName, String chromosome, byte strand, long geneStart, long geneEnd,
            String entrezIds, String karyotypeBand, String synonyms)
    {
        GeneId = geneId;
        GeneName = geneName;
        Chromosome = chromosome;
        Strand = strand;
        GeneStart = geneStart;
        GeneEnd = geneEnd;
        KaryotypeBand = karyotypeBand;

        Synonyms = Arrays.stream(synonyms.split(";")).collect(Collectors.toList());

        if(!entrezIds.equals("NULL"))
        {
            String[] entrezIdStr = entrezIds.split(";");
            EntrezIds = Arrays.stream(entrezIdStr).filter(x -> !x.isEmpty()).map(x -> Integer.parseInt(x)).collect(Collectors.toList());
        }
        else
        {
            EntrezIds = Lists.newArrayList();
        }
    }
}

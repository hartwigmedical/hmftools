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

    public static int GENE_PHASING_REGION_5P_UTR = 0;
    public static int GENE_PHASING_REGION_CODING_0 = 1;
    public static int GENE_PHASING_REGION_CODING_1 = 2;
    public static int GENE_PHASING_REGION_CODING_2 = 3;
    public static int GENE_PHASING_REGION_PROMOTOR = 4;
    public static int GENE_PHASING_REGION_MAX = 5;

    private int mListIndex; // in the set per chromosome
    private int mReverseListIndex;

    private int[] mRegionTotals;

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

        if(!entrezIds.equalsIgnoreCase("NULL"))
        {
            String[] entrezIdStr = entrezIds.split(";");
            EntrezIds = Arrays.stream(entrezIdStr).filter(x -> !x.isEmpty()).map(x -> Integer.parseInt(x)).collect(Collectors.toList());
        }
        else
        {
            EntrezIds = Lists.newArrayList();
        }

        mListIndex= -1;
        mReverseListIndex = -1;
        mRegionTotals = new int[GENE_PHASING_REGION_MAX];
    }

    public int getListIndex() { return mListIndex; }
    public void setListIndex(int index) { mListIndex = index; }

    public int getReverseListIndex() { return mReverseListIndex; }
    public void setReverseListIndex(int index) { mReverseListIndex = index; }

    public int[] getRegionTotals () { return mRegionTotals; }

    public static int phaseToRegion(int phase)
    {
        switch(phase)
        {
            case -1: return GENE_PHASING_REGION_5P_UTR;
            case 0: return GENE_PHASING_REGION_CODING_0;
            case 1: return GENE_PHASING_REGION_CODING_1;
            case 2: return GENE_PHASING_REGION_CODING_2;
        }

        return GENE_PHASING_REGION_5P_UTR;
    }

    public static int regionToPhase(int region)
    {
        if(region == GENE_PHASING_REGION_5P_UTR) return -1;
        if(region == GENE_PHASING_REGION_CODING_0) return 0;
        if(region == GENE_PHASING_REGION_CODING_1) return 1;
        if(region == GENE_PHASING_REGION_CODING_2) return 2;
        if(region == GENE_PHASING_REGION_PROMOTOR) return -1;

        return -1;
    }


}

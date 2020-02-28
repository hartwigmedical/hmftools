package com.hartwig.hmftools.common.variant.structural.annotation;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class EnsemblGeneData
{
    public final String GeneId; // aka StableId
    public final String GeneName;
    public final String Chromosome;
    public final byte Strand;
    public final long GeneStart;
    public final long GeneEnd;
    public final String KaryotypeBand;

    private int mListIndex; // in the set per chromosome
    private int mReverseListIndex;

    public EnsemblGeneData(
            String geneId, String geneName, String chromosome, byte strand, long geneStart, long geneEnd, String karyotypeBand)
    {
        GeneId = geneId;
        GeneName = geneName;
        Chromosome = chromosome;
        Strand = strand;
        GeneStart = geneStart;
        GeneEnd = geneEnd;
        KaryotypeBand = karyotypeBand;

        mListIndex= -1;
        mReverseListIndex = -1;
    }

    public boolean forwardStrand() { return Strand == 1; }
    public boolean reverseStrand() { return Strand == -1; }

    public int getListIndex() { return mListIndex; }
    public void setListIndex(int index) { mListIndex = index; }

    public int getReverseListIndex() { return mReverseListIndex; }
    public void setReverseListIndex(int index) { mReverseListIndex = index; }

}

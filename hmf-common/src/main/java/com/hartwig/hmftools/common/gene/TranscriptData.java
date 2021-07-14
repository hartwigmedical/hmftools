package com.hartwig.hmftools.common.gene;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;

import java.util.List;

import com.google.common.collect.Lists;

public class TranscriptData
{
    public final int TransId;
    public final String TransName;
    public final String GeneId;
    public final boolean IsCanonical;
    public final byte Strand;
    public final int TransStart;
    public final int TransEnd;
    public final Integer CodingStart;
    public final Integer CodingEnd;
    public final String BioType;

    private List<ExonData> mExons;

    public TranscriptData(final int transId, final String transName, final String geneId, final boolean isCanonical, final byte strand,
            int transStart, int transEnd, Integer codingStart, Integer codingEnd, String bioType)
    {
        TransId = transId;
        TransName = transName;
        GeneId = geneId;
        IsCanonical = isCanonical;
        Strand = strand;
        TransStart = transStart;
        TransEnd = transEnd;
        CodingStart = codingStart;
        CodingEnd = codingEnd;
        BioType = bioType;
        mExons = Lists.newArrayList();
    }

    public void setExons(final List<ExonData> exons) { mExons = exons; }
    public List<ExonData> exons() { return mExons; }

    public int length() { return TransEnd - TransStart; }
    public boolean posStrand() { return Strand == POS_STRAND; }

    public String toString()
    {
        return String.format("%d:%s pos(%d-%d) exons(%d) %s",
            TransId, TransName, TransStart, TransEnd, mExons.size(), IsCanonical ? "canonical" : "");
    }
}

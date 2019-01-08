package com.hartwig.hmftools.common.variant.structural.annotation;

public class TranscriptExonData
{
    public final String GeneName;
    public final String TransName;
    public final int TransId;
    public final boolean IsCanonical;
    public final byte Strand;
    public final long TransStart;
    public final long TransEnd;
    public final long ExonStart;
    public final long ExonEnd;
    public final int ExonRank;
    public final int ExonPhase;
    public final int ExonPhaseEnd;
    public final Long CodingStart;
    public final Long CodingEnd;
    public final String BioType;

    public TranscriptExonData(final String geneName, final String transName, final int transId, final boolean isCanonical, final byte strand,
            long transStart, long transEnd, long exonStart, long exonEnd, int exonRank, int exonPhase, int exonPhaseEnd, Long codingStart, Long codingEnd,
            String bioType)
    {
        GeneName = geneName;
        TransName = transName;
        TransId = transId;
        IsCanonical = isCanonical;
        Strand = strand;
        TransStart = transStart;
        TransEnd = transEnd;
        ExonStart = exonStart;
        ExonEnd = exonEnd;
        ExonRank = exonRank;
        ExonPhase = exonPhase;
        ExonPhaseEnd = exonPhaseEnd;
        CodingStart = codingStart;
        CodingEnd = codingEnd;
        BioType = bioType;
    }

    public boolean sameExon(final TranscriptExonData other)
    {
        // checks if an exon is equal regardless of the transcript id
        if(!GeneName.equals(other.GeneName) || Strand != other.Strand)
            return false;

        if(ExonStart != other.ExonStart || ExonEnd != other.ExonEnd)
            return false;

        if(ExonPhase != other.ExonPhase || ExonPhaseEnd != other.ExonPhaseEnd)
            return false;

        return true;
    }

    public boolean overlaps(final TranscriptExonData other)
    {
        // assumes not the same exon
        return !(ExonStart > other.ExonEnd || ExonEnd < other.ExonStart);
    }

}

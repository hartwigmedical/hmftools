package com.hartwig.hmftools.linx.annotators;

public class PseudoGeneMatch
{
    public final String Gene;
    public final int TransId;
    public final String TransName;
    public final int ExonRank;
    public final int ExonLength;
    public int StartHomologyOffset;
    public int EndHomologyOffset;
    public int StartPositionMismatch;
    public int EndPositionMismatch;

    public PseudoGeneMatch(final String gene, final int transId, final String transName, int exonRank, int exonLength)
    {
        Gene = gene;
        TransId = transId;
        TransName = transName;
        ExonRank = exonRank;
        ExonLength = exonLength;
        StartHomologyOffset = 0;
        EndHomologyOffset = 0;
        StartPositionMismatch = 0;
        EndPositionMismatch = 0;
    }

    public boolean isHomologyMatch(boolean isStart)
    {
        return isStart ? StartPositionMismatch == 0 : EndPositionMismatch == 0;
    }

}

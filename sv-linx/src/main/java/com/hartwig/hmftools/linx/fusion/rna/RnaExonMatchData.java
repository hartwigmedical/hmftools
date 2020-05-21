package com.hartwig.hmftools.linx.fusion.rna;

public class RnaExonMatchData
{
    public final int TransId;
    public boolean ExonFound;
    public int ExonRank;
    public int ExonPhase;
    public boolean BoundaryMatch;

    public RnaExonMatchData(final int transId)
    {
        TransId = transId;
        ExonPhase = -1;
        ExonRank = 0;
        ExonFound = false;
        BoundaryMatch = false;
    }

}

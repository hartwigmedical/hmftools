package com.hartwig.hmftools.linx.fusion.rna;

public class RnaExonMatchData
{
    public final String TransName;
    public boolean ExonFound;
    public int ExonRank;
    public int ExonPhase;
    public boolean BoundaryMatch;

    public RnaExonMatchData(final String transName)
    {
        TransName = transName;
        ExonPhase = -1;
        ExonRank = 0;
        ExonFound = false;
        BoundaryMatch = false;
    }

}

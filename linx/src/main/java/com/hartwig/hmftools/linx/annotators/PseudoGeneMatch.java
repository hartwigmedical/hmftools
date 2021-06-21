package com.hartwig.hmftools.linx.annotators;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.seIndex;

public class PseudoGeneMatch
{
    public final String Gene;
    public final int TransId;
    public final String TransName;
    public final int ExonRank;
    public final int ExonLength;
    public int[] HomologyOffset;
    public int[] PositionMismatch;

    public PseudoGeneMatch(final String gene, final int transId, final String transName, int exonRank, int exonLength)
    {
        Gene = gene;
        TransId = transId;
        TransName = transName;
        ExonRank = exonRank;
        ExonLength = exonLength;
        HomologyOffset = new int[SE_PAIR];
        PositionMismatch = new int[SE_PAIR];
    }

    public boolean isHomologyMatch(boolean isStart)
    {
        return PositionMismatch[seIndex(isStart)] == 0;
    }

}

package com.hartwig.hmftools.common.gene;

public class CodingBaseData
{
    public int TotalCodingBases;
    public int CodingBases; // for a given breakend
    public int Phase;

    public static final int PHASE_NONE = -1;
    public static final int PHASE_0 = 0;
    public static final int PHASE_1 = 1;
    public static final int PHASE_2 = 2;

    public CodingBaseData()
    {
        TotalCodingBases = 0;
        CodingBases = 0;
        Phase = PHASE_NONE;
    }

    public CodingBaseData(final int totalCodingBases, final int codingBases, final int phase)
    {
        TotalCodingBases = totalCodingBases;
        CodingBases = codingBases;
        Phase = phase;
    }

}

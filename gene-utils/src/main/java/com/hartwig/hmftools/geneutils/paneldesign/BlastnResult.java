package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

public class BlastnResult
{
    public final String Sequence;
    public final double SumBitScore;
    public final int ResultCount;

    public static final double INVALID_SCORE = -1;
    public static final BlastnResult INVALID_RESULT = new BlastnResult("", INVALID_SCORE, 0);

    public BlastnResult(final String sequence, final double sumBitScore, final int resultCount)
    {
        Sequence = sequence;
        SumBitScore = sumBitScore;
        ResultCount = resultCount;
    }

    public boolean isValid() { return ResultCount > 0 && SumBitScore != INVALID_SCORE; }

    public String toString() { return format("seqLength(%d) score(%.3f) results(%d)", Sequence.length(), SumBitScore, ResultCount); }
}

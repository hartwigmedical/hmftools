package com.hartwig.hmftools.bamtools.biomodalcollapse;

import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.FORWARD_HAIRPIN;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.LOW_QUAL_CUTOFF;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.BiomodalConstants.MISSING_BASE;

import java.util.List;

import org.jetbrains.annotations.Nullable;

public class ReverseComplementMatchInfo
{
    public static final float TOTAL_MISMATCH_PROPORTION_THRESHOLD = 0.3f;
    public static final float HIGH_QUAL_MISMATCH_PROPORTION_THRESHOLD = 0.065f;

    public final int Read1Shift;
    public final int HighQualMismatchCount;
    public final float HighQualMismatchProportion;
    public final int TotalMismatchCount;
    public final float TotalMismatchProportion;

    public ReverseComplementMatchInfo(int read1Shift, int highQualMismatchCount, float highQualMismatchProportion, int totalMismatchCount,
            float totalMismatchProportion)
    {
        Read1Shift = read1Shift;
        HighQualMismatchCount = highQualMismatchCount;
        HighQualMismatchProportion = highQualMismatchProportion;
        TotalMismatchCount = totalMismatchCount;
        TotalMismatchProportion = totalMismatchProportion;
    }

    @Nullable
    public static ReverseComplementMatchInfo findBestReverseComplementMatch(final List<BaseQualPair> seq1,
            final List<BaseQualPair> seq2RevComplement)
    {
        ReverseComplementMatchInfo bestMatch = null;

        int minShift = -seq1.size() + FORWARD_HAIRPIN.length() + 1;
        int maxShift = seq1.size() % 2 == 0 ? seq1.size() / 2 - 1 : seq1.size() / 2;
        for(int read1Shift = minShift; read1Shift <= maxShift; read1Shift++)
        {
            int i1 = 0;
            int i2 = 0;
            if(read1Shift < 0)
            {
                i1 = -read1Shift;
            }

            if(read1Shift > 0)
            {
                i2 = read1Shift;
            }

            int matchCount = 0;
            int lowQualMismatchCount = 0;
            int highQualMismatchCount = 0;
            while(i1 < seq1.size() && i2 < seq2RevComplement.size())
            {
                byte base1 = seq1.get(i1).Base;
                byte base2 = seq2RevComplement.get(i2).Base;
                int qual = min(seq1.get(i1).Qual, seq2RevComplement.get(i2).Qual);
                if(base1 == MISSING_BASE || base2 == MISSING_BASE)
                {
                }
                else if(base1 == base2)
                {
                    matchCount++;
                }
                else if(qual > LOW_QUAL_CUTOFF)
                {
                    highQualMismatchCount++;
                }
                else
                {
                    lowQualMismatchCount++;
                }

                i1++;
                i2++;
            }

            int totalMismatchCount = lowQualMismatchCount + highQualMismatchCount;
            int totalCount = matchCount + totalMismatchCount;

            if(totalCount == 0)
            {
                continue;
            }

            float totalMismatchProp = 1.0f * totalMismatchCount / totalCount;
            float highQualMismatchProp = 1.0f * highQualMismatchCount / totalCount;
            if(totalMismatchProp >= TOTAL_MISMATCH_PROPORTION_THRESHOLD)
            {
                continue;
            }

            if(highQualMismatchProp >= HIGH_QUAL_MISMATCH_PROPORTION_THRESHOLD)
            {
                continue;
            }

            if(bestMatch == null || highQualMismatchProp < bestMatch.HighQualMismatchProportion)
            {
                bestMatch =
                        new ReverseComplementMatchInfo(read1Shift, highQualMismatchCount, highQualMismatchProp, totalMismatchCount, totalMismatchProp);
            }
        }

        return bestMatch;
    }
}

package com.hartwig.hmftools.fastqtools.biomodalcollapse;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.codon.Nucleotides.baseIndex;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalBamUtils.LOW_QUAL_CUTOFF;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalCollapseUtil.getConsensusBase;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalCollapseUtil.modCConsensusBaseQualPair;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalConstants.MISMATCH_BASE;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalConstants.MISSING_BASE;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalConstants.MODC_BASE;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;

public class AlignmentStats
{
    public final int MissingCount;
    public final int MismatchCount;
    public final int HighQualMismatchCountGG;
    public final int HighQualMismatchCountOther;
    public final int LowQualUnambiguousCount;
    public final int LowQualAmbiguousCount;
    public final int ModCGCount;
    public final int ModCOtherCount;
    public final int IndelCount;

    public AlignmentStats(int missingCount, int mismatchCount, int highQualMismatchCountGG, int highQualMismatchCountOther,
            int lowQualUnambiguousCount, int lowQualAmbiguousCount, int modCGCount, int modCOtherCount, int indelCount)
    {
        MissingCount = missingCount;
        MismatchCount = mismatchCount;
        HighQualMismatchCountGG = highQualMismatchCountGG;
        HighQualMismatchCountOther = highQualMismatchCountOther;
        LowQualUnambiguousCount = lowQualUnambiguousCount;
        LowQualAmbiguousCount = lowQualAmbiguousCount;
        ModCGCount = modCGCount;
        ModCOtherCount = modCOtherCount;
        IndelCount = indelCount;
    }

    public static AlignmentStats getAlignmentStats(List<Pair<BaseQualPair, BaseQualPair>> alignment)
    {
        int missingCount = 0;
        int mismatchCount = 0;
        int highQualMismatchCountGG = 0;
        int highQualMismatchCountOther = 0;
        int lowQualUnambiguousCount = 0;
        int lowQualAmbiguousCount = 0;
        int modCGCount = 0;
        int modCOtherCount = 0;
        int indelCount = 0;
        for(int i = 0; i < alignment.size(); i++)
        {
            BaseQualPair base1 = alignment.get(i).getLeft();
            BaseQualPair base2 = alignment.get(i).getRight();
            if(base1 == null || base2 == null)
            {
                indelCount++;
                continue;
            }
            else if(base1.Base == MISSING_BASE || base2.Base == MISSING_BASE)
            {
                missingCount++;
                continue;
            }

            byte consensusBase = base2.Base == (byte) 'G' ? MISMATCH_BASE : getConsensusBase(base1.Base, base2.Base);
            if(baseIndex(consensusBase) != -1)
            {
                continue;
            }

            if(consensusBase == MISMATCH_BASE)
            {
                mismatchCount++;
                if(min(base1.Qual, base2.Qual) > LOW_QUAL_CUTOFF)
                {
                    if(base1.Base == (byte) 'G' && base2.Base == (byte) 'G')
                    {
                        highQualMismatchCountGG++;
                    }
                    else
                    {
                        highQualMismatchCountOther++;
                    }
                }
                else if(base1.Qual > LOW_QUAL_CUTOFF && base1.Base != (byte) 'T')
                {
                    lowQualUnambiguousCount++;
                }
                else if(base2.Qual > LOW_QUAL_CUTOFF && base2.Base != (byte) 'A')
                {
                    lowQualUnambiguousCount++;
                }
                else
                {
                    lowQualAmbiguousCount++;
                }
            }
            else if(consensusBase == MODC_BASE)
            {
                if(i == alignment.size() - 1)
                {
                    modCOtherCount++;
                    continue;
                }

                BaseQualPair nextBase1 = alignment.get(i + 1).getLeft();
                BaseQualPair nextBase2 = alignment.get(i + 1).getRight();
                if(nextBase1 == null || nextBase2 == null)
                {
                    modCOtherCount++;
                    continue;
                }

                if(modCConsensusBaseQualPair(nextBase1, nextBase2).Base == (byte) 'G')
                {
                    modCGCount++;
                    continue;
                }

                modCOtherCount++;
            }
            else
            {
                throw new RuntimeException("Unreachable");
            }
        }

        return new AlignmentStats(missingCount, mismatchCount, highQualMismatchCountGG, highQualMismatchCountOther, lowQualUnambiguousCount, lowQualAmbiguousCount, modCGCount, modCOtherCount, indelCount);
    }
}
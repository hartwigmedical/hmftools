package com.hartwig.hmftools.esvee.assembly.alignment;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.BWA_PENALTY_ADJUST;

import com.hartwig.hmftools.common.bwa.BwaMemAlignParams;
import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.bwa.BwaMemAlignerConfig;

public class AlignerUtils
{
    public static BwaMemAligner createAligner(final String refGenomeImageFile)
    {
        return new BwaMemAligner(createAlignerConfig(refGenomeImageFile));
    }

    private static BwaMemAlignParams createAlignParams()
    {
        BwaMemAlignParams base = BwaMemAlignParams.DEFAULT;
        return base
                .withMismatchPenalty(base.mismatchPenalty() + BWA_PENALTY_ADJUST)
                .withGapOpenPenalty(base.gapOpenPenalty() + BWA_PENALTY_ADJUST);
    }

    private static BwaMemAlignerConfig createAlignerConfig(final String refGenomeImageFile)
    {
        return new BwaMemAlignerConfig(refGenomeImageFile, createAlignParams());
    }
}

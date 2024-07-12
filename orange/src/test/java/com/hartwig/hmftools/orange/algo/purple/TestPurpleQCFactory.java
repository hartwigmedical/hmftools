package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;

import org.jetbrains.annotations.NotNull;

public final class TestPurpleQCFactory
{
    @NotNull
    public static ImmutablePurpleQC.Builder builder()
    {
        return ImmutablePurpleQC.builder()
                .amberMeanDepth(0)
                .unsupportedCopyNumberSegments(0)
                .totalCopyNumberSegments(100)
                .deletedGenes(0)
                .contamination(0D);
    }
}

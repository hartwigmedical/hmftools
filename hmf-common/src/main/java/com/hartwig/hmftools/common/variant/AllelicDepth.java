package com.hartwig.hmftools.common.variant;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;

public interface AllelicDepth {

    int totalReadCount();

    int alleleReadCount();

    default double alleleFrequency() {
        return (double) alleleReadCount() / totalReadCount();
    }

    @NotNull
    static AllelicDepth fromGenotype(@NotNull final Genotype genotype) {
        Preconditions.checkArgument(genotype.hasAD());

        int[] adFields = genotype.getAD();
        final int alleleReadCount = adFields[1];
        int totalReadCount = 0;
        for (final int afField : adFields) {
            totalReadCount += afField;
        }

        return ImmutableAllelicDepthImpl.builder().alleleReadCount(alleleReadCount).totalReadCount(totalReadCount).build();
    }

}

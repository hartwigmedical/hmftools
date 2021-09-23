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

    static boolean containsAllelicDepth(final Genotype genotype) {
        return genotype != null && genotype.hasAD() && genotype.getAD().length > 1;
    }

    @NotNull
    static AllelicDepth fromGenotype(@NotNull final Genotype genotype) {
        Preconditions.checkArgument(genotype.hasAD());
        int[] adFields = genotype.getAD();
        final int alleleReadCount = adFields[1];
        int totalReadCount = totalReadCount(genotype);
        return ImmutableAllelicDepthImpl.builder().alleleReadCount(alleleReadCount).totalReadCount(totalReadCount).build();
    }

    static int totalReadCount(@NotNull final Genotype genotype) {
        // Note: this is a workaround of strelka's DP being only Tier 1
        return genotype.hasDP() ? Math.max(genotype.getDP(), sumReadCount(genotype.getAD())) : sumReadCount(genotype.getAD());
    }

    static int sumReadCount(int[] adFields) {
        int totalReadCount = 0;
        for (final int afField : adFields) {
            totalReadCount += afField;
        }
        return totalReadCount;
    }
}

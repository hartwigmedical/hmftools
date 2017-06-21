package com.hartwig.hmftools.common.gene;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

import org.jetbrains.annotations.NotNull;

class GeneCopyNumberBuilder {

    private final HmfGenomeRegion gene;
    private final ImmutableGeneCopyNumber.Builder builder;

    private double minCopyNumber = Double.MAX_VALUE;
    private double maxCopyNumber = Double.MIN_VALUE;
    private int count;
    private double cumulativeCopyNumber;

    GeneCopyNumberBuilder(@NotNull final HmfGenomeRegion gene) {
        this.gene = gene;
        builder = ImmutableGeneCopyNumber.builder().from(gene);
    }

    void addCopyNumber(@NotNull final PurpleCopyNumber copyNumber) {
        count++;
        minCopyNumber = Math.min(minCopyNumber, copyNumber.averageTumorCopyNumber());
        maxCopyNumber = Math.max(maxCopyNumber, copyNumber.averageTumorCopyNumber());
        cumulativeCopyNumber = gene.overlappingBases(copyNumber) * copyNumber.averageTumorCopyNumber();
    }

    @NotNull
    public String chromosome() {
        return gene.chromosome();
    }

    @NotNull
    public GeneCopyNumber build() {
        return builder.maxCopyNumber(maxCopyNumber)
                .minCopyNumber(minCopyNumber)
                .meanCopyNumber(cumulativeCopyNumber / gene.bases())
                .regions(count)
                .build();
    }
}

package com.hartwig.hmftools.common.gene;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfExonRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

import org.jetbrains.annotations.NotNull;

class GeneCopyNumberBuilder {

    private final ImmutableGeneCopyNumber.Builder builder;

    private double minCopyNumber = Double.MAX_VALUE;
    private double maxCopyNumber = -Double.MAX_VALUE;
    private int somaticCount;
    private int homCount;
    private int het2HomCount;
    private double cumulativeCopyNumber;

    private double previousCopyNumber = -Double.MAX_VALUE;
    private long totalBases;
    private HmfExonRegion exon;
    private PurpleCopyNumber copyNumber;

    GeneCopyNumberBuilder(@NotNull final HmfGenomeRegion gene) {
        builder = ImmutableGeneCopyNumber.builder().from(gene);
    }

    void addExon(@NotNull final HmfExonRegion exon) {
        totalBases += exon.bases();
        this.exon = exon;
        if (copyNumber != null) {
            addOverlap(this.exon, copyNumber);
        }
    }

    void addCopyNumber(@NotNull final PurpleCopyNumber copyNumber) {
        this.copyNumber = copyNumber;
        if (exon != null) {
            addOverlap(exon, this.copyNumber);
        }
    }

    private void addOverlap(@NotNull final HmfExonRegion exon, @NotNull final PurpleCopyNumber copyNumber) {
        long overlap = exon.overlappingBases(copyNumber);
        if (overlap > 0) {
            double currentCopyNumber = copyNumber.averageTumorCopyNumber();

            minCopyNumber = Math.min(minCopyNumber, currentCopyNumber);
            maxCopyNumber = Math.max(maxCopyNumber, currentCopyNumber);
            cumulativeCopyNumber += overlap * currentCopyNumber;

            if (!Doubles.equal(currentCopyNumber, previousCopyNumber)) {
                switch (copyNumber.method()) {
                    case GERMLINE_HOM_DELETION:
                        homCount++;
                        break;
                    case GERMLINE_HET2HOM_DELETION:
                        het2HomCount++;
                        break;
                    default:
                        somaticCount++;
                }
            }

            previousCopyNumber = currentCopyNumber;
        }
    }

    @NotNull
    public GeneCopyNumber build() {
        return builder.maxCopyNumber(maxCopyNumber)
                .minCopyNumber(minCopyNumber)
                .meanCopyNumber(cumulativeCopyNumber / totalBases)
                .somaticRegions(somaticCount)
                .germlineHomRegions(homCount)
                .germlineHet2HomRegions(het2HomCount)
                .build();
    }
}

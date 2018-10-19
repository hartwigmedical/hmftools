package com.hartwig.hmftools.common.purple.qc;

import java.util.List;

import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public final class PurpleQCFactory {

    private PurpleQCFactory() {
    }

    @NotNull
    public static PurpleQC create(@NotNull FittedPurity purity, @NotNull List<PurpleCopyNumber> copyNumbers, @NotNull Gender amberGender,
            @NotNull Gender cobaltGender, @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        int unsupportedSegments = (int) copyNumbers.stream().filter(x -> x.segmentStartSupport() == SegmentSupport.NONE).count();
        int deletedGenes = (int) geneCopyNumbers.stream()
                .filter(x -> !HumanChromosome.fromString(x.chromosome()).equals(HumanChromosome._Y) && x.germlineHet2HomRegions() == 0
                        && x.germlineHomRegions() == 0 && Doubles.lessThan(x.minCopyNumber(), 0.5))
                .count();

        return ImmutablePurpleQC.builder()
                .cobaltGender(cobaltGender)
                .amberGender(amberGender)
                .ploidy(purity.ploidy())
                .unsupportedSegments(unsupportedSegments)
                .deletedGenes(deletedGenes)
                .build();
    }
}

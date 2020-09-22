package com.hartwig.hmftools.common.purple.qc;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.INF;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public final class PurpleQCFactory {

    private static final EnumSet<SegmentSupport> INTERNAL = EnumSet.of(NONE, CENTROMERE, TELOMERE, INF);

    private PurpleQCFactory() {
    }

    @NotNull
    public static PurpleQC create(@NotNull FittedPurity purity, @NotNull List<PurpleCopyNumber> copyNumbers, @NotNull Gender amberGender,
            @NotNull Gender cobaltGender, @NotNull List<GeneCopyNumber> geneCopyNumbers, @NotNull final Set<GermlineAberration> aberrations) {
        boolean containsAnySvSupport = copyNumbers.stream().anyMatch(PurpleQCFactory::isStructuralVariantBreak);

        int unsupportedSegments = containsAnySvSupport ? (int) copyNumbers.stream()
                .filter(x -> x.segmentStartSupport() == NONE && x.segmentEndSupport() == NONE)
                .count() : 0;
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
                .germlineAberrations(aberrations)
                .build();
    }

    private static boolean isStructuralVariantBreak(@NotNull PurpleCopyNumber copyNumber) {
        return !INTERNAL.contains(copyNumber.segmentStartSupport());
    }
}

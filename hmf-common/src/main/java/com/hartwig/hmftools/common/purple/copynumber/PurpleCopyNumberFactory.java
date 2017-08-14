package com.hartwig.hmftools.common.purple.copynumber;

import static java.util.stream.Collectors.toList;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class PurpleCopyNumberFactory {

    @NotNull
    private final PurityAdjuster purityAdjuster;

    private final List<PurpleCopyNumber> highConfidenceRegions;
    private final List<PurpleCopyNumber> smoothedRegions;

    public PurpleCopyNumberFactory(@NotNull final PurityAdjuster purityAdjuster, final List<FittedRegion> fittedRegions) {
        this.purityAdjuster = purityAdjuster;
        smoothedRegions = Lists.newArrayList();
        highConfidenceRegions = Lists.newArrayList();

        final Set<String> orderedChromosomes =
                fittedRegions.stream().map(GenomeRegion::chromosome).collect(Collectors.toCollection(LinkedHashSet::new));

        for (String chromosome : orderedChromosomes) {
            final List<FittedRegion> chromosomeFittedRegions =
                    fittedRegions.stream().filter(matchesChromosome(chromosome)).collect(toList());

            final List<PurpleCopyNumber> highConfidence = highConfidence(chromosomeFittedRegions);
            highConfidenceRegions.addAll(highConfidence);

            final List<FittedRegion> smoothFittedRegions = highConfidence.isEmpty()
                    ? new LowConfidenceSmoothedRegions(purityAdjuster, chromosomeFittedRegions).smoothedRegions()
                    : new HighConfidenceSmoothedRegions(purityAdjuster, highConfidence, chromosomeFittedRegions).smoothedRegions();

            final List<PurpleCopyNumber> copyNumbers = smoothFittedRegions.stream().map(this::create).collect(toList());
            smoothedRegions.addAll(RegionStepFilter.filter(copyNumbers));
        }
    }

    public List<PurpleCopyNumber> highConfidenceRegions() {
        return highConfidenceRegions;
    }

    public List<PurpleCopyNumber> smoothedRegions() {
        return smoothedRegions;
    }

    private List<PurpleCopyNumber> highConfidence(@NotNull final List<FittedRegion> fittedRegions) {
        return new HighConfidenceRegions(purityAdjuster).highConfidence(fittedRegions).stream().map(this::create).collect(toList());
    }

    @NotNull
    private PurpleCopyNumber create(@NotNull final FittedRegion region) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(region.chromosome())
                .start(region.start())
                .end(region.end())
                .bafCount(region.bafCount())
                .averageObservedBAF(region.observedBAF())
                .averageActualBAF(purityAdjuster.purityAdjustedBAF(region.chromosome(), region.tumorCopyNumber(), region.observedBAF()))
                .averageTumorCopyNumber(region.tumorCopyNumber())
                .ratioSupport(region.ratioSupport())
                .structuralVariantSupport(region.structuralVariantSupport())
                .build();
    }

    @NotNull
    private static <T extends GenomeRegion> Predicate<T> matchesChromosome(@NotNull final String chromosome) {
        return t -> t.chromosome().equals(chromosome);
    }
}

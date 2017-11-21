package com.hartwig.hmftools.common.purple.copynumber;

import static java.util.stream.Collectors.toList;

import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.ObservedRegionStatus;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

public class PurpleCopyNumberFactory {

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final List<PurpleCopyNumber> highConfidenceRegions;
    @NotNull
    private final List<PurpleCopyNumber> smoothedRegions;

    public PurpleCopyNumberFactory(boolean experimental, @NotNull final PurityAdjuster purityAdjuster,
            final List<FittedRegion> fittedRegions, final List<StructuralVariant> structuralVariants) {
        this.purityAdjuster = purityAdjuster;
        smoothedRegions = Lists.newArrayList();
        highConfidenceRegions = Lists.newArrayList();

        final Set<String> orderedChromosomes =
                fittedRegions.stream().map(GenomeRegion::chromosome).collect(Collectors.toCollection(LinkedHashSet::new));

        final ListMultimap<String, CombinedRegion> somaticExtentions = ArrayListMultimap.create();

        for (String chromosome : orderedChromosomes) {
            final List<FittedRegion> chromosomeFittedRegions =
                    fittedRegions.stream().filter(matchesChromosome(chromosome)).collect(toList());

            final List<PurpleCopyNumber> highConfidence = highConfidence(chromosomeFittedRegions);
            highConfidenceRegions.addAll(highConfidence);

            final List<FittedRegion> smoothFittedRegions = highConfidence.isEmpty()
                    ? new LowConfidenceSmoothedRegions(purityAdjuster, chromosomeFittedRegions).smoothedRegions()
                    : new HighConfidenceSmoothedRegions(purityAdjuster, highConfidence, chromosomeFittedRegions).smoothedRegions();

            somaticExtentions.putAll(chromosome, SomaticExtension.combinedRegions(purityAdjuster, chromosomeFittedRegions));

            // Old Method
            if (!experimental) {
                smoothedRegions.addAll(RegionStepFilter.filter(toCopyNumber(smoothFittedRegions)));
            }
        }

        if (experimental) {

            final StructuralVariantImplied svImpliedFactory = new StructuralVariantImplied(purityAdjuster);
            final ListMultimap<String, CombinedRegion> svImplied =
                    svImpliedFactory.svImpliedCopyNumber(structuralVariants, somaticExtentions);

            for (HumanChromosome chromosome : HumanChromosome.values()) {
                if (svImplied.containsKey(chromosome.toString())) {
                    smoothedRegions.addAll(toCopyNumber(svImplied.get(chromosome.toString())
                            .stream()
                            .map(CombinedRegion::region)
                            .collect(toList())));
                }
            }
        }

        Collections.sort(smoothedRegions);
    }

    public List<PurpleCopyNumber> highConfidenceRegions() {
        return highConfidenceRegions;
    }

    public List<PurpleCopyNumber> smoothedRegions() {
        return smoothedRegions;
    }

    private List<PurpleCopyNumber> highConfidence(@NotNull final List<FittedRegion> fittedRegions) {
        return new HighConfidenceRegions(purityAdjuster).highConfidence(fittedRegions)
                .stream()
                .map(x -> create(x, SegmentSupport.NONE))
                .collect(toList());
    }

    @NotNull
    private List<PurpleCopyNumber> toCopyNumber(@NotNull final List<FittedRegion> regions) {
        final List<PurpleCopyNumber> result = Lists.newArrayList();
        for (int i = 0; i < regions.size() - 1; i++) {
            final FittedRegion region = regions.get(i);
            final FittedRegion next = regions.get(i + 1);
            result.add(create(region, next.support()));
        }

        if (!regions.isEmpty()) {
            result.add(create(regions.get(regions.size() - 1), SegmentSupport.TELOMERE));
        }

        return PurpleCopyNumberSmoothing.smooth(result);
    }

    @NotNull
    private PurpleCopyNumber create(@NotNull final FittedRegion region, @NotNull final SegmentSupport trailingSupport) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(region.chromosome())
                .start(region.start())
                .end(region.end())
                .bafCount(region.bafCount())
                .averageObservedBAF(region.observedBAF())
                .averageActualBAF(region.tumorBAF())
                .averageTumorCopyNumber(region.tumorCopyNumber())
                .inferred(region.status() == ObservedRegionStatus.CLUSTER)
                .segmentStartSupport(region.support())
                .segmentEndSupport(trailingSupport)
                .build();
    }

    @NotNull
    private static <T extends GenomeRegion> Predicate<T> matchesChromosome(@NotNull final String chromosome) {
        return t -> t.chromosome().equals(chromosome);
    }
}

package com.hartwig.hmftools.common.purple.copynumber;

import static java.util.stream.Collectors.toList;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class PurpleCopyNumberFactory {

    private final double purity;

    private final List<PurpleCopyNumber> highConfidenceRegions;
    private final List<PurpleCopyNumber> smoothedRegions;

    public PurpleCopyNumberFactory(final double purity, final List<FittedRegion> fittedRegions) {
        this.purity = purity;
        smoothedRegions = Lists.newArrayList();
        highConfidenceRegions = Lists.newArrayList();

        final Set<String> orderedChromosomes = fittedRegions.stream().map(GenomeRegion::chromosome).collect(Collectors.toCollection(LinkedHashSet::new));

        for (String chromosome : orderedChromosomes) {
            final List<FittedRegion> chromosomeFittedRegions = fittedRegions.stream().filter(matchesChromosome(chromosome)).collect(toList());

            final List<PurpleCopyNumber> highConfidence = highConfidence(chromosome, chromosomeFittedRegions);
            highConfidenceRegions.addAll(highConfidence);

            final List<PurpleCopyNumber> smooth = new SmoothedRegions(purity, highConfidence, chromosomeFittedRegions).getSmoothedRegions();
            smoothedRegions.addAll(RegionStepFilter.filter(smooth));
        }
    }

    public List<PurpleCopyNumber> highConfidenceRegions() {
        return highConfidenceRegions;
    }

    public List<PurpleCopyNumber> smoothedRegions() {
        return smoothedRegions;
    }

    private List<PurpleCopyNumber> highConfidence(final String chromosome, @NotNull final List<FittedRegion> fittedRegions) {
        return new HighConfidenceRegions(purity).highConfidence(fittedRegions);
    }


    @NotNull
    private static <T extends GenomeRegion> Predicate<T> matchesChromosome(@NotNull final String chromosome) {
        return t -> t.chromosome().equals(chromosome);
    }
}

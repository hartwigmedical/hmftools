package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.Optional;
import java.util.function.Function;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantLegCopyNumberFactory<T extends GenomeRegion> {

    @NotNull
    private final Function<T, Double> copyNumberExtractor;

    public StructuralVariantLegCopyNumberFactory(@NotNull final Function<T, Double> copyNumberExtractor) {
        this.copyNumberExtractor = copyNumberExtractor;
    }

    @NotNull
    public StructuralVariantLegCopyNumber create(@NotNull final StructuralVariantLeg leg,
            @NotNull final Multimap<Chromosome, T> copyNumbers) {
        return create(leg, GenomeRegionSelectorFactory.createImproved(copyNumbers));

    }

    @NotNull
    public StructuralVariantLegCopyNumber create(@NotNull final StructuralVariantLeg leg, @NotNull final GenomeRegionSelector<T> selector) {
        final GenomePosition svPositionLeft = GenomePositions.create(leg.chromosome(), leg.cnaPosition() - 1);
        final GenomePosition svPositionRight = GenomePositions.create(leg.chromosome(), leg.cnaPosition());
        final Optional<Double> left =
                selector.select(svPositionLeft).flatMap(x -> Optional.ofNullable(copyNumberExtractor.apply(x))).map(x -> Math.max(0, x));
        final Optional<Double> right =
                selector.select(svPositionRight).flatMap(x -> Optional.ofNullable(copyNumberExtractor.apply(x))).map(x -> Math.max(0, x));

        return ImmutableStructuralVariantLegCopyNumberImpl.builder()
                .from(leg)
                .leftCopyNumber(left)
                .rightCopyNumber(right)
                .build();
    }
}

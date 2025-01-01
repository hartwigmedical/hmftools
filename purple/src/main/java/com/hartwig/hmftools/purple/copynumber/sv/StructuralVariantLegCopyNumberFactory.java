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
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

import org.jetbrains.annotations.NotNull;

public class StructuralVariantLegCopyNumberFactory<T extends GenomeRegion>
{
    private final Function<T, Double> mCopyNumberExtractor;

    public StructuralVariantLegCopyNumberFactory(final Function<T, Double> copyNumberExtractor)
    {
        mCopyNumberExtractor = copyNumberExtractor;
    }

    public StructuralVariantLegCopyNumber create(final StructuralVariantLeg leg, final Multimap<Chromosome, T> copyNumbers)
    {
        return create(leg, GenomeRegionSelectorFactory.createImproved(copyNumbers));
    }

    public StructuralVariantLegCopyNumber create(final StructuralVariantLeg leg, final GenomeRegionSelector<T> selector)
    {
        GenomePosition svPositionLeft = GenomePositions.create(leg.chromosome(), leg.cnaPosition() - 1);
        GenomePosition svPositionRight = GenomePositions.create(leg.chromosome(), leg.cnaPosition());

        Optional<Double> left =
                selector.select(svPositionLeft).flatMap(x -> Optional.ofNullable(mCopyNumberExtractor.apply(x))).map(x -> Math.max(0, x));

        Optional<Double> right =
                selector.select(svPositionRight).flatMap(x -> Optional.ofNullable(mCopyNumberExtractor.apply(x))).map(x -> Math.max(0, x));

        return ImmutableStructuralVariantLegCopyNumberImpl.builder()
                .from(leg)
                .leftCopyNumber(left)
                .rightCopyNumber(right)
                .build();
    }
}

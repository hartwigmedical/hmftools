package com.hartwig.hmftools.purple.copynumber.sv;

import java.util.Collections;
import java.util.List;
import java.util.Optional;
import java.util.function.Function;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class SvLegPloidyFactory<T extends GenomeRegion>
{
    private static final double VAF_TO_USE_READ_DEPTH = 0.75;
    private static final double RECIPROCAL_VAF_TO_USE_READ_DEPTH = reciprocalVAF(VAF_TO_USE_READ_DEPTH);

    private final double mAverageCopyNumber;
    private final int mAverageReadDepth;

    private final PurityAdjuster mPurityAdjuster;
    private final SvLegCopyNumberFactory<T> mCopyNumberFactory;

    public SvLegPloidyFactory(final PurityAdjuster purityAdjuster, final Function<T, Double> copyNumberExtractor)
    {
        this(0, 0, purityAdjuster, copyNumberExtractor);
    }

    public SvLegPloidyFactory(int averageReadDepth, double averageCopyNumber, final PurityAdjuster purityAdjuster,
            final Function<T, Double> copyNumberExtractor)
    {
        mAverageCopyNumber = averageCopyNumber;
        mAverageReadDepth = averageReadDepth;
        mPurityAdjuster = purityAdjuster;
        mCopyNumberFactory = new SvLegCopyNumberFactory<>(copyNumberExtractor);
    }

    public List<StructuralVariantLegPloidy> create(final StructuralVariant variant, final Multimap<Chromosome, T> copyNumbers)
    {
        final List<StructuralVariantLegPloidy> result = Lists.newArrayList();
        final List<StructuralVariantLegs> allLegs = SvLegsFactory.create(variant);

        for(StructuralVariantLegs leg : allLegs)
        {
            result.addAll(create(leg, copyNumbers));
        }

        Collections.sort(result);
        return result;
    }

    public List<StructuralVariantLegPloidy> create(final List<StructuralVariant> variants, final Multimap<Chromosome, T> copyNumbers)
    {
        final List<StructuralVariantLegPloidy> result = Lists.newArrayList();
        final List<StructuralVariantLegs> allLegs = SvLegsFactory.create(variants);

        for(StructuralVariantLegs leg : allLegs)
        {
            result.addAll(create(leg, copyNumbers));
        }

        Collections.sort(result);
        return result;
    }

    public List<StructuralVariantLegPloidy> create(final StructuralVariantLegs legs, final Multimap<Chromosome, T> copyNumbers)
    {
        Optional<StructuralVariantLegPloidy> start = legs.start() != null ?
                create(legs.start(), GenomeRegionSelectorFactory.createImproved(copyNumbers)) : Optional.empty();

        Optional<StructuralVariantLegPloidy> end = legs.end() != null ?
                create(legs.end(), GenomeRegionSelectorFactory.createImproved(copyNumbers)) : Optional.empty();

        if(!start.isPresent() && !end.isPresent())
        {
            return Collections.emptyList();
        }

        List<StructuralVariantLegPloidy> result = Lists.newArrayList();
        double startWeight = start.map(StructuralVariantLegPloidy::weight).orElse(0D);
        double startPloidy = start.map(StructuralVariantLegPloidy::unweightedImpliedPloidy).orElse(0D);
        double endWeight = end.map(StructuralVariantLegPloidy::weight).orElse(0D);
        double endPloidy = end.map(StructuralVariantLegPloidy::unweightedImpliedPloidy).orElse(0D);

        double totalWeight = startWeight + endWeight;
        double averagePloidy = (startWeight * startPloidy + endWeight * endPloidy) / totalWeight;

        if(start.isPresent())
        {
            start.get().setWeight(totalWeight);
            start.get().setAverageImpliedPloidy(averagePloidy);
            result.add(start.get());
        }

        if(end.isPresent())
        {
            end.get().setWeight(totalWeight);
            end.get().setAverageImpliedPloidy(averagePloidy);
            result.add(end.get());
        }

        Collections.sort(result);
        return result;
    }

    public Optional<StructuralVariantLegPloidy> create(final StructuralVariantLeg leg,
            final GenomeRegionSelector<T> selector)
    {
        final StructuralVariantLegCopyNumber legCopyNumber = mCopyNumberFactory.create(leg, selector);
        return create(leg, legCopyNumber.leftCopyNumber(), legCopyNumber.rightCopyNumber());
    }

    private Optional<StructuralVariantLegPloidy> create(
            final StructuralVariantLeg leg, final Double leftCopyNumber, final Double rightCopyNumber)
    {
        final Double largerCopyNumber;
        final Double smallerCopyNumber;

        if(leg.orientation() == 1)
        {
            largerCopyNumber = leftCopyNumber;
            smallerCopyNumber = rightCopyNumber;
        }
        else
        {
            largerCopyNumber = rightCopyNumber;
            smallerCopyNumber = leftCopyNumber;
        }

        if(largerCopyNumber== null && smallerCopyNumber == null)
        {
            return Optional.empty();
        }

        final Double observedVaf = leg.alleleFrequency();
        if(observedVaf == null)
        {
            return Optional.empty();
        }

        final double adjustedVaf;
        final double ploidy;
        final double weight;

        if(largerCopyNumber != null)
        {
            double copyNumber = largerCopyNumber;
            adjustedVaf = mPurityAdjuster.purityAdjustedVAF(leg.chromosome(), Math.max(0.001, copyNumber), observedVaf);
            ploidy = adjustedVaf * copyNumber;
            weight = 1;
        }
        else
        {
            double copyNumber = smallerCopyNumber;
            double reciprocalVAF = reciprocalVAF(observedVaf);
            adjustedVaf = mPurityAdjuster.purityAdjustedVAF(leg.chromosome(), Math.max(0.001, copyNumber), reciprocalVAF);

            if(mAverageReadDepth > 0 && (Double.isInfinite(adjustedVaf) || Doubles.greaterThan(adjustedVaf,
                    RECIPROCAL_VAF_TO_USE_READ_DEPTH)))
            {
                final Integer tumorVariantFragmentCount = leg.tumorVariantFragmentCount();
                if(tumorVariantFragmentCount != null && tumorVariantFragmentCount > 0)
                {
                    ploidy = readDepthImpliedPloidy(tumorVariantFragmentCount);
                }
                else
                {
                    return Optional.empty();
                }

            }
            else if(!Double.isFinite(reciprocalVAF))
            {
                return Optional.empty();
            }
            else
            {
                ploidy = adjustedVaf * copyNumber;
            }

            weight = 1 / (1 + Math.pow(Math.max(copyNumber, 2) / Math.min(Math.max(copyNumber, 0.01), 2), 2));
        }

        StructuralVariantLegCopyNumber svLeg = new StructuralVariantLegCopyNumber(leg, leftCopyNumber, rightCopyNumber);

        StructuralVariantLegPloidy legPloidy = new StructuralVariantLegPloidy(
                svLeg, observedVaf, adjustedVaf, 0, ploidy, weight);

        return Optional.of(legPloidy);
    }

    private static double reciprocalVAF(double vaf)
    {
        return vaf / (1 - vaf);
    }

    private double readDepthImpliedPloidy(int tumorVariantFragmentCount)
    {
        return mAverageCopyNumber * tumorVariantFragmentCount / mAverageReadDepth;
    }
}

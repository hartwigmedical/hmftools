package com.hartwig.hmftools.purple.somatic;

import java.util.List;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.purple.region.FittedRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.variant.AllelicDepth;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class PurityAdjustedSomaticVariantFactory
{
    private final PurityAdjuster mPurityAdjuster;
    private final GenomeRegionSelector<PurpleCopyNumber> mCopyNumberSelector;
    private final GenomeRegionSelector<FittedRegion> mFittedRegionSelector;
    private final String mSample;

    public PurityAdjustedSomaticVariantFactory(
            final String sample, final PurityAdjuster purityAdjuster,
            final List<PurpleCopyNumber> copyNumbers, final List<FittedRegion> fittedRegions)
    {
        this(sample, purityAdjuster, Multimaps.fromRegions(copyNumbers), Multimaps.fromRegions(fittedRegions));
    }

    private PurityAdjustedSomaticVariantFactory(
            final String sample, final PurityAdjuster purityAdjuster,
            final Multimap<Chromosome, PurpleCopyNumber> copyNumbers, final Multimap<Chromosome, FittedRegion> fittedRegions)
    {
        mSample = sample;
        mPurityAdjuster = purityAdjuster;
        mCopyNumberSelector = GenomeRegionSelectorFactory.createImproved(copyNumbers);
        mFittedRegionSelector = GenomeRegionSelectorFactory.createImproved(fittedRegions);
    }

    public VariantContext enrich(final VariantContext variant)
    {
        processVariant(variant);
        return variant;
    }

    public void processVariant(final VariantContext variant)
    {
        final Genotype genotype = variant.getGenotype(mSample);
        if(genotype != null && genotype.hasAD() && HumanChromosome.contains(variant.getContig()))
        {
            final GenomePosition position = GenomePositions.create(variant.getContig(), variant.getStart());
            final AllelicDepth depth = AllelicDepth.fromGenotype(genotype);
            final PurityAdjustedSomaticVariantBuilder builder = PurityAdjustedSomaticVariantBuilder.fromVariantContext(variant);
            enrich(position, depth, builder);
        }
    }

    private void enrich(final GenomePosition position, final AllelicDepth depth, final PurityAdjustedSomaticVariantBuilder builder)
    {
        mCopyNumberSelector.select(position).ifPresent(x -> applyPurityAdjustment(x, depth, builder));
        mFittedRegionSelector.select(position).ifPresent(x -> builder.germlineStatus(x.germlineStatus()));
    }

    private void applyPurityAdjustment(
            final PurpleCopyNumber purpleCopyNumber, final AllelicDepth depth, final PurityAdjustedSomaticVariantBuilder builder)
    {
        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
        double vaf = mPurityAdjuster.purityAdjustedVAF(purpleCopyNumber.chromosome(), Math.max(0.001, copyNumber), depth.alleleFrequency());
        double ploidy = Math.max(0, vaf * copyNumber);

        boolean biallelic = Doubles.lessOrEqual(copyNumber, 0) || Doubles.greaterOrEqual(ploidy, copyNumber - 0.5);

        builder.adjustedCopyNumber(copyNumber)
                .adjustedVAF(vaf)
                .variantCopyNumber(ploidy)
                .biallelic(biallelic)
                .minorAlleleCopyNumber(purpleCopyNumber.minorAlleleCopyNumber());
    }
}

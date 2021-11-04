package com.hartwig.hmftools.common.variant;

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
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Multimaps;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class PurityAdjustedSomaticVariantFactory {

    @NotNull
    private final PurityAdjuster purityAdjuster;
    @NotNull
    private final GenomeRegionSelector<PurpleCopyNumber> copyNumberSelector;
    @NotNull
    private final GenomeRegionSelector<FittedRegion> fittedRegionSelector;
    @NotNull
    private final String sample;

    public PurityAdjustedSomaticVariantFactory(@NotNull String sample, @NotNull PurityAdjuster purityAdjuster,
            @NotNull final List<PurpleCopyNumber> copyNumbers, @NotNull final List<FittedRegion> fittedRegions) {
        this(sample, purityAdjuster, Multimaps.fromRegions(copyNumbers), Multimaps.fromRegions(fittedRegions));
    }

    private PurityAdjustedSomaticVariantFactory(@NotNull String sample, @NotNull PurityAdjuster purityAdjuster,
            @NotNull final Multimap<Chromosome, PurpleCopyNumber> copyNumbers,
            @NotNull final Multimap<Chromosome, FittedRegion> fittedRegions) {
        this.sample = sample;
        this.purityAdjuster = purityAdjuster;
        this.copyNumberSelector = GenomeRegionSelectorFactory.createImproved(copyNumbers);
        this.fittedRegionSelector = GenomeRegionSelectorFactory.createImproved(fittedRegions);
    }

    @NotNull
    public VariantContext enrich(@NotNull final VariantContext variant) {
        final Genotype genotype = variant.getGenotype(sample);
        if (genotype != null && genotype.hasAD() && HumanChromosome.contains(variant.getContig())) {
            final GenomePosition position = GenomePositions.create(variant.getContig(), variant.getStart());
            final AllelicDepth depth = AllelicDepth.fromGenotype(genotype);
            enrich(position, depth, PurityAdjustedSomaticVariantBuilder.fromVariantContext(variant));
        }
        return variant;
    }

    private void enrich(@NotNull final GenomePosition position, @NotNull final AllelicDepth depth,
            @NotNull final PurityAdjustedSomaticVariantBuilder builder) {
        copyNumberSelector.select(position).ifPresent(x -> applyPurityAdjustment(x, depth, builder));
        fittedRegionSelector.select(position).ifPresent(x -> builder.germlineStatus(x.germlineStatus()));
    }

    private void applyPurityAdjustment(@NotNull final PurpleCopyNumber purpleCopyNumber, @NotNull final AllelicDepth depth,
            @NotNull final PurityAdjustedSomaticVariantBuilder builder) {
        double copyNumber = purpleCopyNumber.averageTumorCopyNumber();
        double vaf = purityAdjuster.purityAdjustedVAF(purpleCopyNumber.chromosome(), Math.max(0.001, copyNumber), depth.alleleFrequency());
        double ploidy = Math.max(0, vaf * copyNumber);

        boolean biallelic = Doubles.lessOrEqual(copyNumber, 0) || Doubles.greaterOrEqual(ploidy, copyNumber - 0.5);

        builder.adjustedCopyNumber(copyNumber)
                .adjustedVAF(vaf)
                .variantCopyNumber(ploidy)
                .biallelic(biallelic)
                .minorAlleleCopyNumber(purpleCopyNumber.minorAlleleCopyNumber());
    }
}

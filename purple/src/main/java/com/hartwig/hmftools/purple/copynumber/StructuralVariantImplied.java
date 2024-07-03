package com.hartwig.hmftools.purple.copynumber;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidyFactory;
import com.hartwig.hmftools.common.sv.StructuralVariant;

public class StructuralVariantImplied
{
    private final StructuralVariantLegPloidyFactory<CombinedRegion> mSvPloidyFactory;

    public StructuralVariantImplied(int averageReadDepth, double averageCopyNumber, final PurityAdjuster purityAdjuster)
    {
        mSvPloidyFactory = new StructuralVariantLegPloidyFactory<>(averageReadDepth,
                averageCopyNumber,
                purityAdjuster,
                x -> x.isProcessed() ? x.tumorCopyNumber() : 0);
    }

    public ListMultimap<Chromosome, CombinedRegion> svImpliedCopyNumber(
            final List<StructuralVariant> structuralVariants, final ListMultimap<Chromosome, CombinedRegion> copyNumbers)
    {
        long previousMissingCopyNumbers = copyNumbers.size();
        long currentMissingCopyNumbers = missingCopyNumberCount(copyNumbers);

        while(currentMissingCopyNumbers < previousMissingCopyNumbers && currentMissingCopyNumbers > 0)
        {
            final ListMultimap<Chromosome, CombinedRegion> processedCopyNumbers = ArrayListMultimap.create();
            for(Chromosome chromosome : copyNumbers.keySet())
            {
                copyNumbers.get(chromosome)
                        .stream()
                        .filter(CombinedRegion::isProcessed)
                        .forEach(x -> processedCopyNumbers.put(chromosome, x));
            }

            final List<StructuralVariantLegPloidy> ploidyList = mSvPloidyFactory.create(structuralVariants, processedCopyNumbers);
            final Map<GenomePosition, StructuralVariantLegPloidy> ploidyMap =
                    ploidyList.stream().collect(Collectors.toMap(x -> GenomePositions.create(x.chromosome(), x.cnaPosition()), x -> x));

            for(Chromosome chromosome : HumanChromosome.values())
            {
                final List<CombinedRegion> chromosomeCopyNumbers = copyNumbers.get(chromosome);
                boolean svInferred = false;
                for(final CombinedRegion copyNumber : chromosomeCopyNumbers)
                {
                    if(implyCopyNumberFromSV(copyNumber))
                    {
                        final Optional<StructuralVariantLegPloidy> optionalStart =
                                Optional.ofNullable(ploidyMap.get(GenomePositions.create(copyNumber.chromosome(), copyNumber.start())))
                                        .filter(x -> x.impliedRightCopyNumberWeight() > 0);
                        final Optional<StructuralVariantLegPloidy> optionalEnd =
                                Optional.ofNullable(ploidyMap.get(GenomePositions.create(copyNumber.chromosome(), copyNumber.end() + 1)))
                                        .filter(x -> x.impliedLeftCopyNumberWeight() > 0);
                        if(optionalStart.isPresent() || optionalEnd.isPresent())
                        {
                            svInferred = true;
                            inferCopyNumberFromStructuralVariants(copyNumber, optionalStart, optionalEnd);
                        }
                    }
                }

                // Extend structural variant segments
                if(svInferred)
                {
                    ExtendStructuralVariant.extendStructuralVariants(chromosomeCopyNumbers);
                }
            }

            previousMissingCopyNumbers = currentMissingCopyNumbers;
            currentMissingCopyNumbers = missingCopyNumberCount(copyNumbers);
        }

        return copyNumbers;
    }

    private void inferCopyNumberFromStructuralVariants(final CombinedRegion region,
            final Optional<StructuralVariantLegPloidy> start, final Optional<StructuralVariantLegPloidy> end)
    {
        region.resetDepthWindowCount();
        region.setTumorCopyNumber(CopyNumberMethod.STRUCTURAL_VARIANT, inferCopyNumberFromStructuralVariants(start, end));
    }

    @VisibleForTesting
    static double inferCopyNumberFromStructuralVariants(
            final Optional<StructuralVariantLegPloidy> start, final Optional<StructuralVariantLegPloidy> end)
    {
        final double startWeight = start.map(StructuralVariantLegPloidy::impliedRightCopyNumberWeight).orElse(0d);
        final double startCopyNumber = start.map(StructuralVariantLegPloidy::impliedRightCopyNumber).orElse(0d);

        final double endWeight = end.map(StructuralVariantLegPloidy::impliedLeftCopyNumberWeight).orElse(0d);
        final double endCopyNumber = end.map(StructuralVariantLegPloidy::impliedLeftCopyNumber).orElse(0d);

        double unconstrainedResult = (startCopyNumber * startWeight + endCopyNumber * endWeight) / (startWeight + endWeight);
        return Math.max(0, unconstrainedResult);
    }

    private long missingCopyNumberCount(Multimap<?, CombinedRegion> copyNumbers)
    {
        return copyNumbers.values().stream().filter(this::implyCopyNumberFromSV).count();
    }

    private boolean implyCopyNumberFromSV(final CombinedRegion copyNumber)
    {
        return !copyNumber.isProcessed();
    }
}

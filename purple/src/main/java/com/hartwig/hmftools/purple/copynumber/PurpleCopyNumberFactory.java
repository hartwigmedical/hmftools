package com.hartwig.hmftools.purple.copynumber;

import static java.util.stream.Collectors.toList;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.purple.purity.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.jetbrains.annotations.NotNull;

public class PurpleCopyNumberFactory
{
    private final List<PurpleCopyNumber> mSomaticCopyNumbers = Lists.newArrayList();

    private final double mPloidy;
    private final int mAverageReadDepth;
    private final int mMinTumorRatioCount;
    private final int mMinTumorRatioCountAtCentromere;
    private final PurityAdjuster mPurityAdjuster;
    private final CobaltChromosomes mCobaltChromosomes;

    public PurpleCopyNumberFactory(
            int minTumorRatioCount, int minTumorRatioCountAtCentromere, int averageReadDepth, double ploidy,
            @NotNull final PurityAdjuster purityAdjuster, @NotNull final CobaltChromosomes cobaltChromosomes)
    {
        mPurityAdjuster = purityAdjuster;
        mMinTumorRatioCount = minTumorRatioCount;
        mMinTumorRatioCountAtCentromere = minTumorRatioCountAtCentromere;
        mAverageReadDepth = averageReadDepth;
        mPloidy = ploidy;
        mCobaltChromosomes = cobaltChromosomes;
    }

    public void invoke(final List<ObservedRegion> fittedRegions, final List<StructuralVariant> structuralVariants)
    {
        mSomaticCopyNumbers.clear();

        final ExtendDiploid extendDiploid =
                new ExtendDiploid(new AlleleTolerance(mPurityAdjuster), mMinTumorRatioCount, mMinTumorRatioCountAtCentromere);

        final PopulateUnknown populateUnknownFactory = new PopulateUnknown(mCobaltChromosomes);

        final ListMultimap<Chromosome, CombinedRegion> diploidExtension = ArrayListMultimap.create();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            final List<ObservedRegion> chromosomeFittedRegions =
                    fittedRegions.stream().filter(matchesChromosome(chromosome)).collect(toList());

            final List<CombinedRegion> diploidExtended = extendDiploid.extendDiploid(chromosomeFittedRegions);
            final List<CombinedRegion> nonDiploidExtended = ExtendNonDiploid.nonDiploid(diploidExtended);

            diploidExtension.putAll(chromosome, nonDiploidExtended);
        }

        final StructuralVariantImplied svImpliedFactory = new StructuralVariantImplied(mAverageReadDepth, mPloidy, mPurityAdjuster);
        final ListMultimap<Chromosome, CombinedRegion> allSVImplied =
                svImpliedFactory.svImpliedCopyNumber(structuralVariants, diploidExtension);

        for(final HumanChromosome chromosome : HumanChromosome.values())
        {
            final ExtendDiploidBAF extendDiploidBAF = new ExtendDiploidBAF(simpleVariants(chromosome, structuralVariants));

            final List<CombinedRegion> svImplied = allSVImplied.get(chromosome);
            final List<CombinedRegion> longArmExtended = ExtendLongArm.extendLongArm(svImplied);
            final List<CombinedRegion> populateUnknown = populateUnknownFactory.populateUnknown(longArmExtended);
            final List<CombinedRegion> somatics = extendDiploidBAF.extendBAF(populateUnknown);

            mSomaticCopyNumbers.addAll(toCopyNumber(somatics));
        }
    }

    @NotNull
    private List<StructuralVariant> simpleVariants(HumanChromosome chromosome, final List<StructuralVariant> structuralVariants)
    {
        return structuralVariants.stream().filter(x ->
        {
            StructuralVariantLeg end = x.end();
            return end != null && HumanChromosome.contains(x.start().chromosome()) && HumanChromosome.contains(end.chromosome())
                    && x.start().chromosome().equals(end.chromosome()) && HumanChromosome.fromString(x.start().chromosome())
                    .equals(chromosome);
        }).collect(toList());
    }

    @NotNull
    public List<PurpleCopyNumber> copyNumbers()
    {
        return mSomaticCopyNumbers;
    }

    @NotNull
    private static List<PurpleCopyNumber> toCopyNumber(final List<CombinedRegion> regions)
    {
        final List<PurpleCopyNumber> result = Lists.newArrayList();
        for(int i = 0; i < regions.size() - 1; i++)
        {
            final CombinedRegion region = regions.get(i);
            final CombinedRegion next = regions.get(i + 1);
            result.add(toCopyNumber(region, next.region().support()));
        }

        if(!regions.isEmpty())
        {
            result.add(toCopyNumber(regions.get(regions.size() - 1), SegmentSupport.TELOMERE));
        }

        return result;
    }

    @NotNull
    private static PurpleCopyNumber toCopyNumber(final CombinedRegion region, final SegmentSupport trailingSupport)
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(region.chromosome())
                .start(region.start())
                .end(region.end())
                .bafCount(region.bafCount())
                .averageObservedBAF(region.region().observedBAF())
                .averageActualBAF(region.tumorBAF())
                .averageTumorCopyNumber(region.tumorCopyNumber())
                .segmentStartSupport(region.region().support())
                .segmentEndSupport(trailingSupport)
                .method(region.copyNumberMethod())
                .depthWindowCount(region.region().depthWindowCount())
                .gcContent(region.region().gcContent())
                .minStart(region.region().minStart())
                .maxStart(region.region().maxStart())
                .build();
    }

    @NotNull
    private static <T extends GenomeRegion> Predicate<T> matchesChromosome(final Chromosome chromosome)
    {
        return t -> HumanChromosome.fromString(t.chromosome()).equals(chromosome);
    }
}

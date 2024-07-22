package com.hartwig.hmftools.purple.copynumber;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.CDKN2A_DELETION_REGION;
import static com.hartwig.hmftools.purple.copynumber.ExtendUtils.populateUnknown;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.purple.region.ObservedRegion;

public class PurpleCopyNumberFactory
{
    private final List<PurpleCopyNumber> mSomaticCopyNumbers;

    private final double mPloidy;
    private final int mAverageReadDepth;
    private final int mMinTumorRatioCount;
    private final int mMinTumorRatioCountAtCentromere;
    private final PurityAdjuster mPurityAdjuster;
    private final CobaltChromosomes mCobaltChromosomes;

    public PurpleCopyNumberFactory(
            int minTumorRatioCount, int minTumorRatioCountAtCentromere, int averageReadDepth, double ploidy,
            final PurityAdjuster purityAdjuster, final CobaltChromosomes cobaltChromosomes)
    {
        mPurityAdjuster = purityAdjuster;
        mMinTumorRatioCount = minTumorRatioCount;
        mMinTumorRatioCountAtCentromere = minTumorRatioCountAtCentromere;
        mAverageReadDepth = averageReadDepth;
        mPloidy = ploidy;
        mCobaltChromosomes = cobaltChromosomes;
        mSomaticCopyNumbers = Lists.newArrayList();
    }

    public void buildCopyNumbers(final List<ObservedRegion> fittedRegions, final List<StructuralVariant> structuralVariants)
    {
        mSomaticCopyNumbers.clear();

        ExtendDiploid extendDiploid = new ExtendDiploid(new AlleleTolerance(mPurityAdjuster), mMinTumorRatioCount, mMinTumorRatioCountAtCentromere);

        ListMultimap<Chromosome, CombinedRegion> diploidExtension = ArrayListMultimap.create();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            List<ObservedRegion> chromosomeFittedRegions = fittedRegions.stream()
                    .filter(x -> chromosome.matches(x.chromosome())).collect(toList());

            List<CombinedRegion> diploidExtended = extendDiploid.extendDiploid(chromosomeFittedRegions);
            List<CombinedRegion> nonDiploidExtended = ExtendNonDiploid.nonDiploid(diploidExtended);

            diploidExtension.putAll(chromosome, nonDiploidExtended);
        }

        StructuralVariantImplied svImpliedFactory = new StructuralVariantImplied(mAverageReadDepth, mPloidy, mPurityAdjuster);

        ListMultimap<Chromosome, CombinedRegion> allSVImplied = svImpliedFactory.svImpliedCopyNumber(structuralVariants, diploidExtension);

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            ExtendDiploidBAF extendDiploidBAF = new ExtendDiploidBAF(simpleVariants(chromosome, structuralVariants));

            List<CombinedRegion> svImplied = allSVImplied.get(chromosome);
            List<CombinedRegion> longArmExtended = ExtendLongArm.extendLongArm(svImplied);
            List<CombinedRegion> populateUnknown = populateUnknown(longArmExtended, mCobaltChromosomes);
            List<CombinedRegion> somatics = extendDiploidBAF.extendBAF(populateUnknown);

            mSomaticCopyNumbers.addAll(toCopyNumber(somatics));
        }
    }

    private List<StructuralVariant> simpleVariants(final HumanChromosome chromosome, final List<StructuralVariant> structuralVariants)
    {
        return structuralVariants.stream().filter(x ->
        {
            StructuralVariantLeg end = x.end();
            return end != null && HumanChromosome.contains(x.start().chromosome()) && HumanChromosome.contains(end.chromosome())
                    && x.start().chromosome().equals(end.chromosome()) && HumanChromosome.fromString(x.start().chromosome())
                    .equals(chromosome);
        }).collect(toList());
    }

    public List<PurpleCopyNumber> copyNumbers() { return mSomaticCopyNumbers; }

    private static List<PurpleCopyNumber> toCopyNumber(final List<CombinedRegion> regions)
    {
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();

        for(int i = 0; i < regions.size(); i++)
        {
            final CombinedRegion region = regions.get(i);

            int copyNumberStartPos = region.start();

            if(i > 0)
            {
                ObservedRegion firstRegion = region.regions().get(0);
                PurpleCopyNumber prevCopyNumber = copyNumbers.get(copyNumbers.size() - 1);

                if(firstRegion.germlineStatus() == GermlineStatus.DIPLOID && firstRegion.minStart() < firstRegion.maxStart()
                && firstRegion.minStart() == prevCopyNumber.end() + 1)
                {
                    copyNumberStartPos = (firstRegion.minStart() + firstRegion.maxStart()) / 2;

                    PurpleCopyNumber newPrevCopyNumber = ImmutablePurpleCopyNumber.builder()
                            .from(prevCopyNumber)
                            .end(copyNumberStartPos - 1)
                            .build();

                    copyNumbers.set(copyNumbers.size() - 1, newPrevCopyNumber);
                }
            }

            final SegmentSupport trailingSupport = i < regions.size() - 1 ? regions.get(i + 1).region().support() : SegmentSupport.TELOMERE;

            PurpleCopyNumber copyNumber = toCopyNumber(region, copyNumberStartPos, trailingSupport);
            copyNumbers.add(copyNumber);
        }

        return copyNumbers;
    }

    private static PurpleCopyNumber toCopyNumber(final CombinedRegion region, int copyNumberStartPos, final SegmentSupport trailingSupport)
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(region.chromosome())
                .start(copyNumberStartPos)
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

    public static boolean validateCopyNumbers(final List<PurpleCopyNumber> copyNumbers)
    {
        boolean isValid = true;

        for(int i = 1; i < copyNumbers.size(); ++i)
        {
            PurpleCopyNumber copyNumber = copyNumbers.get(i);

            /*
            if(!positionsWithin(copyNumber.minStart(), copyNumber.maxStart(), copyNumber.start(), copyNumber.end()))
            {
                PPL_LOGGER.error("purple copy-number({}:{}-{}) has invalid min/maxStart({}-{})",
                        copyNumber.chromosome(), copyNumber.start(), copyNumber.end(),
                        copyNumber.minStart(), copyNumber.maxStart());

                isValid = false;
            }
            */

            PurpleCopyNumber prevCopyNumber = copyNumbers.get(i - 1);

            if(!copyNumber.chromosome().equals(prevCopyNumber.chromosome()))
                continue;

            if(copyNumber.start() <= prevCopyNumber.end())
            {
                PPL_LOGGER.error("purple copy-number({}:{}-{}) overlaps previous({}-{})",
                        copyNumber.chromosome(), copyNumber.start(), copyNumber.end(), prevCopyNumber.start(), prevCopyNumber.end());
                isValid = false;
            }
        }

        return isValid;
    }

    private static final HumanChromosome CDKN2A_CHR = HumanChromosome.fromString(CDKN2A_DELETION_REGION.Chromosome);

    public static double calculateDeletedDepthWindows(final List<PurpleCopyNumber> copyNumbers)
    {
        int totalDepthWindows = 0;
        int deletedDepthWindows = 0;

        // excluding chr Y and a large region around CDKN2A
        for(PurpleCopyNumber copyNumber : copyNumbers)
        {
            totalDepthWindows += copyNumber.depthWindowCount();

            HumanChromosome chromosome = HumanChromosome.fromString(copyNumber.chromosome());

            if(chromosome == HumanChromosome._Y)
                continue;

            if(copyNumber.averageTumorCopyNumber() >= 0.5)
                continue;

            int deletedWindows = copyNumber.depthWindowCount();

            if(chromosome == CDKN2A_CHR
            && positionsOverlap(copyNumber.start(), copyNumber.end(), CDKN2A_DELETION_REGION.start(), CDKN2A_DELETION_REGION.end()))
            {
                int baseOverlap = min(copyNumber.end(), CDKN2A_DELETION_REGION.end()) - max(copyNumber.start(), CDKN2A_DELETION_REGION.start());
                double nonGeneDeletedFraction = (copyNumber.length() - baseOverlap) / (double)copyNumber.length();
                deletedWindows = (int)round(nonGeneDeletedFraction * deletedWindows);
            }

            deletedDepthWindows += deletedWindows;
        }

        return totalDepthWindows > 0 ? deletedDepthWindows / (double)totalDepthWindows : 0;
    }
}

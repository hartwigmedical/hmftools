package com.hartwig.hmftools.common.purple.copynumber;

import static java.util.stream.Collectors.toList;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.tolerance.AlleleTolerance;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.jetbrains.annotations.NotNull;

public class PurpleCopyNumberFactory {

    private final List<PurpleCopyNumber> somaticCopyNumbers = Lists.newArrayList();
    private final List<PurpleCopyNumber> germlineDeletions = Lists.newArrayList();

    private final double ploidy;
    private final int averageReadDepth;
    private final Gender gender;
    private final int minTumorRatioCount;
    private final int minTumorRatioCountAtCentromere;
    private final PurityAdjuster purityAdjuster;

    public PurpleCopyNumberFactory(Gender gender, int minTumorRatioCount, int minTumorRatioCountAtCentromere, int averageReadDepth, double ploidy,
            @NotNull final PurityAdjuster purityAdjuster) {
        this.purityAdjuster = purityAdjuster;
        this.gender = gender;
        this.minTumorRatioCount = minTumorRatioCount;
        this.minTumorRatioCountAtCentromere = minTumorRatioCountAtCentromere;
        this.averageReadDepth = averageReadDepth;
        this.ploidy = ploidy;

    }

    public void invoke(final List<FittedRegion> fittedRegions, final List<StructuralVariant> structuralVariants) {
        somaticCopyNumbers.clear();
        germlineDeletions.clear();

        final ExtractGermlineDeletions extendGermline = new ExtractGermlineDeletions(gender);
        final ExtendDiploid extendDiploid =
                new ExtendDiploid(new AlleleTolerance(purityAdjuster), minTumorRatioCount, minTumorRatioCountAtCentromere);
        final PopulateUnknown populateUnknownFactory = new PopulateUnknown(gender);

        final ListMultimap<Chromosome, CombinedRegion> diploidExtension = ArrayListMultimap.create();
        for (HumanChromosome chromosome : HumanChromosome.values()) {
            final List<FittedRegion> chromosomeFittedRegions = fittedRegions.stream().filter(matchesChromosome(chromosome)).collect(toList());

            final List<CombinedRegion> diploidExtended = extendDiploid.extendDiploid(chromosomeFittedRegions);
            final List<CombinedRegion> nonDiploidExtended = ExtendNonDiploid.nonDiploid(diploidExtended);

            diploidExtension.putAll(chromosome, nonDiploidExtended);
        }

        final StructuralVariantImplied svImpliedFactory = new StructuralVariantImplied(averageReadDepth, ploidy, purityAdjuster);
        final ListMultimap<Chromosome, CombinedRegion> allSVImplied =
                svImpliedFactory.svImpliedCopyNumber(structuralVariants, diploidExtension);

        for (final HumanChromosome chromosome : HumanChromosome.values()) {
            final ExtendDiploidBAF extendDiploidBAF = new ExtendDiploidBAF(simpleVariants(chromosome, structuralVariants));

            final List<CombinedRegion> svImplied = allSVImplied.get(chromosome);
            final List<CombinedRegion> longArmExtended = ExtendLongArm.extendLongArm(svImplied);
            final List<CombinedRegion> populateUnknown = populateUnknownFactory.populateUnknown(longArmExtended);
            final List<CombinedRegion> somatics = extendDiploidBAF.extendBAF(populateUnknown);

            final List<CombinedRegion> germlineDeletions = extendGermline.extractGermlineDeletions(somatics);

            this.somaticCopyNumbers.addAll(toCopyNumber(somatics));
            this.germlineDeletions.addAll(germlineDeletions.stream().map(x -> toCopyNumber(x, SegmentSupport.UNKNOWN)).collect(toList()));
        }
    }

    @NotNull
    private List<StructuralVariant> simpleVariants(HumanChromosome chromosome, final List<StructuralVariant> structuralVariants) {
        return structuralVariants.stream().filter(x -> {
            StructuralVariantLeg end = x.end();
            return end != null && HumanChromosome.contains(x.start().chromosome()) && HumanChromosome.contains(end.chromosome())
                    && x.start().chromosome().equals(end.chromosome()) && HumanChromosome.fromString(x.start().chromosome())
                    .equals(chromosome);
        }).collect(toList());
    }

    @NotNull
    public List<PurpleCopyNumber> copyNumbers() {
        return somaticCopyNumbers;
    }

    @NotNull
    public List<PurpleCopyNumber> germlineDeletions() {
        return germlineDeletions;
    }

    @NotNull
    private static List<PurpleCopyNumber> toCopyNumber(@NotNull final List<CombinedRegion> regions) {
        final List<PurpleCopyNumber> result = Lists.newArrayList();
        for (int i = 0; i < regions.size() - 1; i++) {
            final CombinedRegion region = regions.get(i);
            final CombinedRegion next = regions.get(i + 1);
            result.add(toCopyNumber(region, next.region().support()));
        }

        if (!regions.isEmpty()) {
            result.add(toCopyNumber(regions.get(regions.size() - 1), SegmentSupport.TELOMERE));
        }

        return result;
    }

    @NotNull
    private static PurpleCopyNumber toCopyNumber(@NotNull final CombinedRegion region, @NotNull final SegmentSupport trailingSupport) {
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
    private static <T extends GenomeRegion> Predicate<T> matchesChromosome(@NotNull final Chromosome chromosome) {
        return t -> HumanChromosome.fromString(t.chromosome()).equals(chromosome);
    }
}

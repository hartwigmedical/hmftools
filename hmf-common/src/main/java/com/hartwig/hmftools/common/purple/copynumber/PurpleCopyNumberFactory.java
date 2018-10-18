package com.hartwig.hmftools.common.purple.copynumber;

import static java.util.stream.Collectors.toList;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.tolerance.AlleleTolerance;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.jetbrains.annotations.NotNull;

public class PurpleCopyNumberFactory {

    private final List<PurpleCopyNumber> somaticCopyNumbers = Lists.newArrayList();
    private final List<PurpleCopyNumber> germlineDeletions = Lists.newArrayList();

    private final Gender gender;
    private final int minTumorRatioCount;
    private final int minTumorRatioCountAtCentromere;
    private final PurityAdjuster purityAdjuster;

    public PurpleCopyNumberFactory(int minTumorRatioCount, int minTumorRatioCountAtCentromere,
            @NotNull final PurityAdjuster purityAdjuster) {
        this.purityAdjuster = purityAdjuster;
        this.gender = purityAdjuster.gender();
        this.minTumorRatioCount = minTumorRatioCount;
        this.minTumorRatioCountAtCentromere = minTumorRatioCountAtCentromere;

    }

    public void invoke(final List<FittedRegion> fittedRegions, final List<StructuralVariant> structuralVariants) {
        somaticCopyNumbers.clear();
        germlineDeletions.clear();

        final ExtendGermline extendGermline = new ExtendGermline(gender);
        final ExtendDiploid extendDiploid =
                new ExtendDiploid(new AlleleTolerance(purityAdjuster), minTumorRatioCount, minTumorRatioCountAtCentromere);
        final PopulateUnknown populateUnknownFactory = new PopulateUnknown(gender);

        final ListMultimap<String, CombinedRegion> diploidExtension = ArrayListMultimap.create();
        for (HumanChromosome chromosome : HumanChromosome.values()) {
            final List<FittedRegion> chromosomeFittedRegions =
                    fittedRegions.stream().filter(matchesChromosome(chromosome)).collect(toList());
            diploidExtension.putAll(chromosome.toString(), extendDiploid.extendDiploid(chromosomeFittedRegions));
        }

        final StructuralVariantImplied svImpliedFactory = new StructuralVariantImplied(purityAdjuster);
        final ListMultimap<String, CombinedRegion> allSVImplied =
                svImpliedFactory.svImpliedCopyNumber(structuralVariants, diploidExtension);

        for (HumanChromosome chromosome : HumanChromosome.values()) {
            if (allSVImplied.containsKey(chromosome.toString())) {
                final List<CombinedRegion> svImplied = allSVImplied.get(chromosome.toString());
                final List<CombinedRegion> longArmExtended = ExtendLongArm.extendLongArm(svImplied);
                final List<CombinedRegion> populateUnknown = populateUnknownFactory.populateUnknown(longArmExtended);
                final List<CombinedRegion> bafExtended = ExtendDiploidBAF.extendBAF(populateUnknown);

                final List<CombinedRegion> somatics = extendGermline.extendGermlineAmplifications(bafExtended);
                final List<CombinedRegion> germlineDeletions = extendGermline.extractGermlineDeletions(bafExtended);

                this.somaticCopyNumbers.addAll(toCopyNumber(somatics));
                this.germlineDeletions.addAll(germlineDeletions.stream()
                        .map(x -> toCopyNumber(x, SegmentSupport.UNKNOWN))
                        .collect(toList()));
            }
        }
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

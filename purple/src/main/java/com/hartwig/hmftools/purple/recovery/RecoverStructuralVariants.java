package com.hartwig.hmftools.purple.recovery;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.INFERRED;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.RECOVERED;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.RECOVERY_FILTER;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.RECOVERY_METHOD;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.SVTYPE;

import java.io.Closeable;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegCopyNumberChangeFactory;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidyFactory;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class RecoverStructuralVariants implements Closeable
{
    private static final Allele REF_ALLELE = Allele.create("N", true);
    private static final Allele INCREASING_ALLELE = Allele.create(".N", false);
    private static final Allele DECREASING_ALLELE = Allele.create("N.", false);

    static final int UNBALANCED_MIN_DEPTH_WINDOW_COUNT = 5;
    private static final double UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE = 0.6;
    private static final double UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE_AS_PERCENT_OF_COPY_NUMBER = 0.2;

    private final PurityAdjuster mPurityAdjuster;
    private final ListMultimap<Chromosome, PurpleCopyNumber> mAllCopyNumbers;
    private final StructuralVariantLegPloidyFactory<PurpleCopyNumber> mPloidyFactory;
    private final RecoveredVariantFactory mRecoveredVariantFactory;

    private int mCounter = 0;

    public RecoverStructuralVariants(@NotNull final PurityAdjuster purityAdjuster, @NotNull final String recoveryVCF,
            @NotNull final List<PurpleCopyNumber> allCopyNumbers)
    {
        this(purityAdjuster, new RecoveredVariantFactory(purityAdjuster, recoveryVCF), allCopyNumbers);
    }

    RecoverStructuralVariants(@NotNull final PurityAdjuster purityAdjuster, @NotNull final RecoveredVariantFactory factory,
            @NotNull final List<PurpleCopyNumber> allCopyNumbers)
    {
        mPurityAdjuster = purityAdjuster;
        mAllCopyNumbers = Multimaps.fromRegions(allCopyNumbers);
        mPloidyFactory = new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);
        mRecoveredVariantFactory = factory;
    }

    @NotNull
    public Collection<VariantContext> recoverVariants(@NotNull final List<StructuralVariant> currentVariants) throws IOException
    {
        final Map<String, VariantContext> result = Maps.newHashMap();

        recoverFromUnexplainedSegments().forEach(x -> result.put(x.getID(), x));
        recoverFromUnbalancedVariants(currentVariants, result.values()).forEach(x -> result.put(x.getID(), x));

        return result.values();
    }

    @VisibleForTesting
    @NotNull
    List<VariantContext> recoverFromUnbalancedVariants(
            @NotNull final List<StructuralVariant> currentVariants,
            @NotNull final Collection<VariantContext> recovered) throws IOException
    {
        final StructuralVariantLegCopyNumberChangeFactory changeFactory =
                new StructuralVariantLegCopyNumberChangeFactory(mPurityAdjuster, mAllCopyNumbers, currentVariants);

        final List<VariantContext> result = Lists.newArrayList();
        for(final StructuralVariant variant : currentVariants)
        {
            int recoverCount = 0;
            boolean attemptRecovery = false;
            boolean unbalancedStart = false;
            boolean unbalancedEnd = false;

            final List<StructuralVariantLegPloidy> legs = mPloidyFactory.create(variant, mAllCopyNumbers);
            for(StructuralVariantLegPloidy leg : legs)
            {
                boolean isStart = leg.position() == variant.start().position();

                double copyNumber = leg.adjustedCopyNumber();
                double copyNumberChange = changeFactory.copyNumberChange(leg);
                double expectedCopyNumberChange = leg.averageImpliedPloidy();
                double unexplainedCopyNumberChange = Math.max(0, expectedCopyNumberChange - copyNumberChange);

                if(isUnbalanced(unexplainedCopyNumberChange, copyNumber) && !isCloseToRecoveredVariant(leg, recovered))
                {
                    if(isStart)
                    {
                        unbalancedStart = true;
                    }
                    else
                    {
                        unbalancedEnd = true;
                    }

                    final List<PurpleCopyNumber> chromosomeCopyNumbers = mAllCopyNumbers.get(HumanChromosome.fromString(leg.chromosome()));
                    int index = indexOf(leg.cnaPosition(), chromosomeCopyNumbers);

                    if(index > 0)
                    {
                        final PurpleCopyNumber prev = chromosomeCopyNumbers.get(index - 1);
                        final PurpleCopyNumber current = chromosomeCopyNumbers.get(index);

                        if(current.segmentStartSupport() != SegmentSupport.MULTIPLE && isSupportedByDepthWindowCounts(prev, current))
                        {
                            attemptRecovery = true;
                            int expectedOrientation = -1 * leg.orientation();
                            final Optional<RecoveredVariant> optionalRecoveredVariant = mRecoveredVariantFactory.recoverVariantAtIndex(
                                    expectedOrientation,
                                    unexplainedCopyNumberChange,
                                    index,
                                    chromosomeCopyNumbers);
                            if(optionalRecoveredVariant.isPresent())
                            {
                                final RecoveredVariant recoveredVariant = optionalRecoveredVariant.get();
                                result.addAll(toContext(recoveredVariant, "UNBALANCED_SV"));
                                recoverCount++;
                            }
                        }
                    }
                }
            }

            if(unbalancedStart != unbalancedEnd && attemptRecovery && recoverCount == 0)
            {
                final StructuralVariantLeg leg = unbalancedStart ? variant.start() : variant.end();
                assert (leg != null);
                result.add(infer(leg));
            }
        }

        return result;
    }

    private static boolean isCloseToRecoveredVariant(StructuralVariantLegPloidy legPloidy, Collection<VariantContext> recovered)
    {
        for(VariantContext other : recovered)
        {
            if(legPloidy.chromosome().equals(other.getContig())
            && Math.abs(legPloidy.position() - other.getStart()) <= UNBALANCED_MIN_DEPTH_WINDOW_COUNT * 1000)
            {
                return true;
            }
        }

        return false;
    }

    @NotNull
    private List<VariantContext> recoverFromUnexplainedSegments() throws IOException
    {
        final List<VariantContext> result = Lists.newArrayList();

        for(Chromosome chromosome : mAllCopyNumbers.keySet())
        {
            final List<PurpleCopyNumber> chromosomeCopyNumbers = mAllCopyNumbers.get(chromosome);

            for(int index = 1; index < chromosomeCopyNumbers.size() - 1; index++)
            {
                final PurpleCopyNumber current = chromosomeCopyNumbers.get(index);
                if(current.segmentStartSupport() == SegmentSupport.NONE)
                {
                    PurpleCopyNumber prev = chromosomeCopyNumbers.get(index - 1);
                    double unexplainedCopyNumberChange = Math.abs(prev.averageTumorCopyNumber() - current.averageTumorCopyNumber());

                    int expectedOrientation = Doubles.greaterThan(current.averageTumorCopyNumber(), prev.averageTumorCopyNumber()) ? -1 : 1;

                    final Optional<RecoveredVariant> optionalRecoveredVariant = mRecoveredVariantFactory.recoverVariantAtIndex(
                            expectedOrientation,
                            unexplainedCopyNumberChange,
                            index,
                            chromosomeCopyNumbers);
                    optionalRecoveredVariant.ifPresent(recoveredVariant -> result.addAll(toContext(recoveredVariant,
                            "UNSUPPORTED_BREAKEND")));
                }
            }
        }

        return result;
    }

    private static <T extends GenomeRegion> int indexOf(long cnaPosition, @NotNull final List<T> regions)
    {
        assert (!regions.isEmpty());
        for(int i = 0; i < regions.size(); i++)
        {
            if(regions.get(i).start() == cnaPosition)
            {
                return i;
            }
        }

        return -1;
    }

    @NotNull
    private static List<VariantContext> toContext(@NotNull final RecoveredVariant variant, @NotNull final String partialMethod)
    {
        final List<VariantContext> result = Lists.newArrayList();
        final Set<String> recoveryFilterSet = filterSet(variant.context());
        final VariantContext mate = variant.mate();

        final String recoveryMethod;
        if(mate != null)
        {
            final GenomePosition matePosition = GenomePositions.create(mate.getContig(), mate.getStart());
            final GenomePosition contextPosition = GenomePositions.create(variant.context().getContig(), variant.context().getStart());
            final boolean contextIsStart = contextPosition.compareTo(matePosition) <= 0;
            recoveryFilterSet.addAll(filterSet(mate));

            recoveryMethod = partialMethod + (contextIsStart ? "_START" : "_END");
            result.add(addRecoveryDetails(mate, recoveryMethod, recoveryFilterSet.stream().sorted().collect(Collectors.toList())));

        }
        else
        {
            recoveryMethod = partialMethod + "_START";
        }

        result.add(addRecoveryDetails(variant.context(), recoveryMethod, recoveryFilterSet.stream().sorted().collect(Collectors.toList())));
        return result;
    }

    @NotNull
    private VariantContext infer(@NotNull final StructuralVariantLeg leg)
    {
        // Note: Opposite orientation to leg!
        final Allele allele = leg.orientation() < 0 ? DECREASING_ALLELE : INCREASING_ALLELE;
        final Collection<Allele> alleles = Lists.newArrayList(REF_ALLELE, allele);

        return new VariantContextBuilder("purple", leg.chromosome(), leg.position(), leg.position(), alleles).filter(INFERRED)
                .id("unbalanced_" + mCounter++)
                .attribute(SVTYPE, "BND")
                .attribute(INFERRED, true)
                .attribute(RECOVERED, true)
                .attribute(RECOVERY_METHOD, "UNBALANCED_SV_START")
                .noGenotypes()
                .make();
    }

    @NotNull
    private static Set<String> filterSet(@NotNull VariantContext variantContext)
    {
        return variantContext.isNotFiltered() ? Sets.newHashSet("PASS") : Sets.newHashSet(variantContext.getFilters());
    }

    @NotNull
    private static VariantContext addRecoveryDetails(
            @NotNull final VariantContext context, @NotNull final String recoveryMethod, @NotNull final List<String> recoveryFilters)
    {
        return new VariantContextBuilder(context).unfiltered()
                .attribute(RECOVERED, true)
                .attribute(RECOVERY_METHOD, recoveryMethod)
                .attribute(RECOVERY_FILTER, recoveryFilters)
                .make();
    }

    private static boolean isUnbalanced(double unexplainedCopyNumberChange, double copyNumber)
    {
        return Doubles.greaterOrEqual(unexplainedCopyNumberChange,
                UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE_AS_PERCENT_OF_COPY_NUMBER * copyNumber) && Doubles.greaterOrEqual(
                unexplainedCopyNumberChange,
                UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE);
    }

    private static boolean isSupportedByDepthWindowCounts(@NotNull final PurpleCopyNumber prev, @Nullable final PurpleCopyNumber next)
    {
        return prev.depthWindowCount() >= UNBALANCED_MIN_DEPTH_WINDOW_COUNT && (next == null
                || next.depthWindowCount() >= UNBALANCED_MIN_DEPTH_WINDOW_COUNT);
    }

    @Override
    public void close() throws IOException
    {
        mRecoveredVariantFactory.close();
    }
}

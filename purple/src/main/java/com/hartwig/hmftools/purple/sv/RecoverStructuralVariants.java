package com.hartwig.hmftools.purple.sv;

import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED;
import static com.hartwig.hmftools.common.sv.SvVcfTags.RECOVERED;
import static com.hartwig.hmftools.common.sv.SvVcfTags.RECOVERY_FILTER;
import static com.hartwig.hmftools.common.sv.SvVcfTags.RECOVERY_METHOD;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SVTYPE;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.RECOVERY_UNBALANCED_MIN_DEPTH_WINDOW_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.RECOVERY_UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE;
import static com.hartwig.hmftools.purple.PurpleConstants.RECOVERY_UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE_PERC;

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
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.SampleData;
import com.hartwig.hmftools.purple.SampleDataFiles;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegCopyNumberChangeFactory;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidyFactory;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.utils.collection.Multimaps;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class RecoverStructuralVariants implements Closeable
{
    private static final Allele REF_ALLELE = Allele.create("N", true);
    private static final Allele INCREASING_ALLELE = Allele.create(".N", false);
    private static final Allele DECREASING_ALLELE = Allele.create("N.", false);

    private final PurityAdjuster mPurityAdjuster;
    private final ListMultimap<Chromosome, PurpleCopyNumber> mAllCopyNumbers;
    private final StructuralVariantLegPloidyFactory<PurpleCopyNumber> mPloidyFactory;
    private final RecoveredVariantFactory mRecoveredVariantFactory;

    private int mCounter = 0;

    public RecoverStructuralVariants(
            final PurpleConfig config, final PurityAdjuster purityAdjuster, final String recoveryVCF, final List<PurpleCopyNumber> allCopyNumbers)
    {
        RecoveredVariantFactory svFactory = new RecoveredVariantFactory(purityAdjuster, recoveryVCF, config);

        mPurityAdjuster = purityAdjuster;
        mAllCopyNumbers = Multimaps.fromRegions(allCopyNumbers);
        mPloidyFactory = new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);
        mRecoveredVariantFactory = svFactory;
    }

    RecoverStructuralVariants(
            final PurityAdjuster purityAdjuster, final RecoveredVariantFactory factory, final List<PurpleCopyNumber> allCopyNumbers)
    {
        mPurityAdjuster = purityAdjuster;
        mAllCopyNumbers = Multimaps.fromRegions(allCopyNumbers);
        mPloidyFactory = new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);
        mRecoveredVariantFactory = factory;
    }

    public static int recoverStructuralVariants(
            final SampleData sampleData, final SampleDataFiles sampleDataFiles, final PurpleConfig config,
            final PurityAdjuster purityAdjuster, final List<PurpleCopyNumber> copyNumbers)
    {
        if(sampleDataFiles.RecoveredSvVcfFile.isEmpty())
            return 0;

        PPL_LOGGER.info("loading recovery candidates from {}", sampleDataFiles.RecoveredSvVcfFile);

        RecoverStructuralVariants recovery = new RecoverStructuralVariants(
                config, purityAdjuster, sampleDataFiles.RecoveredSvVcfFile, copyNumbers);

        try
        {
            final Collection<VariantContext> recoveredVariants = recovery.recoverVariants(sampleData.SvCache.variants());

            if(!recoveredVariants.isEmpty())
            {
                recoveredVariants.forEach(x -> sampleData.SvCache.addVariant(x));
            }

            return recoveredVariants.size();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to load recovery SVs: {}", e.toString());
            return 0;
        }
    }

    public Collection<VariantContext> recoverVariants(final List<StructuralVariant> currentVariants) throws IOException
    {
        final Map<String, VariantContext> result = Maps.newHashMap();

        recoverFromUnexplainedSegments().forEach(x -> result.put(x.getID(), x));
        recoverFromUnbalancedVariants(currentVariants, result.values()).forEach(x -> result.put(x.getID(), x));

        return result.values();
    }

    @VisibleForTesting
    List<VariantContext> recoverFromUnbalancedVariants(
            final List<StructuralVariant> currentVariants, final Collection<VariantContext> recovered) throws IOException
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

                            final List<RecoveredVariant> candidates = mRecoveredVariantFactory.recoverVariantAtIndex(
                                    expectedOrientation, unexplainedCopyNumberChange, index, chromosomeCopyNumbers);

                            final Optional<RecoveredVariant> optionalRecoveredVariant = mRecoveredVariantFactory.findTopCandidate(candidates);

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
            && Math.abs(legPloidy.position() - other.getStart()) <= RECOVERY_UNBALANCED_MIN_DEPTH_WINDOW_COUNT * 1000)
            {
                return true;
            }
        }

        return false;
    }

    public List<VariantContext> recoverFromUnexplainedSegments() throws IOException
    {
        final List<VariantContext> result = Lists.newArrayList();

        final List<List<RecoveredVariant>> candidateLists = Lists.newArrayList();

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

                    final List<RecoveredVariant> candidates = mRecoveredVariantFactory.recoverVariantAtIndex(
                            expectedOrientation, unexplainedCopyNumberChange, index, chromosomeCopyNumbers);

                    if(!candidates.isEmpty())
                        candidateLists.add(candidates);
                }
            }
        }

        // now prioritise across candidate locations (ie by CN index
        for(int i = 0; i < candidateLists.size(); ++i)
        {
            final List<RecoveredVariant> variants1 = candidateLists.get(i);

            if(variants1.isEmpty())
                continue;

            // check for other lists with a matching variant
            RecoveredVariant topVariant = null;
            boolean topHasRecoveredMate = false;
            int topMateListIndex = -1;

            for(RecoveredVariant variant : variants1)
            {
                boolean foundRecoveredMate = false;
                int mateListIndex = -1;

                if(variant.mate() != null)
                {
                    for(int j = i + 1; j < candidateLists.size(); ++j)
                    {
                        final List<RecoveredVariant> variants2 = candidateLists.get(j);

                        if(variants2.stream()
                                .filter(x -> x.mate() != null)
                                .anyMatch(x -> x.mate().getID().equals(variant.context().getID())))
                        {
                            foundRecoveredMate = true;
                            mateListIndex = j;
                            break;
                        }
                    }
                }

                if(topVariant == null || foundRecoveredMate && !topHasRecoveredMate
                || variant.context().getPhredScaledQual() > topVariant.context().getPhredScaledQual())
                {
                    topVariant = variant;
                    topHasRecoveredMate = foundRecoveredMate;
                    topMateListIndex = mateListIndex;
                }
            }

            result.addAll(toContext(topVariant, "UNSUPPORTED_BREAKEND"));

            // remove other candidates from any matched mate's list
            if(topMateListIndex >= 0)
                candidateLists.get(topMateListIndex).clear();
        }

        return result;
    }

    private static <T extends GenomeRegion> int indexOf(long cnaPosition, final List<T> regions)
    {
        assert (!regions.isEmpty());
        for(int i = 0; i < regions.size(); i++)
        {
            if(regions.get(i).start() == cnaPosition)
                return i;
        }

        return -1;
    }

    private static List<VariantContext> toContext(final RecoveredVariant variant, final String partialMethod)
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

    private VariantContext infer(final StructuralVariantLeg leg)
    {
        // Note: Opposite orientation to leg!
        final Allele allele = leg.orientation() < 0 ? DECREASING_ALLELE : INCREASING_ALLELE;
        final Collection<Allele> alleles = Lists.newArrayList(REF_ALLELE, allele);

        return new VariantContextBuilder("purple", leg.chromosome(), leg.position(), leg.position(), alleles).filter(INFERRED)
                .id("unbalanced_" + mCounter++)
                .attribute(SVTYPE, StructuralVariantType.BND.toString())
                .attribute(INFERRED, true)
                .attribute(RECOVERED, true)
                .attribute(RECOVERY_METHOD, "UNBALANCED_SV_START")
                .noGenotypes()
                .make();
    }

    private static Set<String> filterSet(VariantContext variantContext)
    {
        return variantContext.isNotFiltered() ? Sets.newHashSet("PASS") : Sets.newHashSet(variantContext.getFilters());
    }

    private static VariantContext addRecoveryDetails(
            final VariantContext context, final String recoveryMethod, final List<String> recoveryFilters)
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
                RECOVERY_UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE_PERC * copyNumber)
                && Doubles.greaterOrEqual(unexplainedCopyNumberChange, RECOVERY_UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE);
    }

    private static boolean isSupportedByDepthWindowCounts(final PurpleCopyNumber prev, @Nullable final PurpleCopyNumber next)
    {
        return prev.depthWindowCount() >= RECOVERY_UNBALANCED_MIN_DEPTH_WINDOW_COUNT
                && (next == null || next.depthWindowCount() >= RECOVERY_UNBALANCED_MIN_DEPTH_WINDOW_COUNT);
    }

    @Override
    public void close() throws IOException
    {
        mRecoveredVariantFactory.close();
    }
}

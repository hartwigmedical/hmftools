package com.hartwig.hmftools.common.variant.recovery;

import static java.util.Comparator.comparingDouble;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.RECOVERY_FILTER;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.RECOVERY_METHOD;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.Closeable;
import java.io.IOException;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.StringJoiner;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.position.GenomePositions;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegCopyNumberChangeFactory;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidyFactory;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.vcf.VCFCodec;

public class RecoverStructuralVariants implements Closeable {

    private static final double MIN_LENGTH = 1000;
    private static final double MIN_MATE_QUAL_SCORE = 350;
    private static final double MIN_SINGLE_QUAL_SCORE = 1000;
    private static final double UNBALANCED_MIN_DEPTH_WINDOW_COUNT = 5;
    private static final double UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE = 0.5;
    private static final double UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE_AS_PERCENT_OF_COPY_NUMBER = 0.2;

    private static final double MIN_PLOIDY = 0.5;
    private static final double MIN_PLOIDY_AS_PERCENTAGE_OF_COPY_NUMBER_CHANGE = 0.5;

    private static final int MIN_MATE_UNCERTAINTY = 150;
    private static final String AF_FILTERED = "af";

    private static final Comparator<RecoveredVariant> QUALITY_COMPARATOR = comparingDouble(x -> x.context().getPhredScaledQual());

    private final PurityAdjuster purityAdjuster;
    private final ListMultimap<Chromosome, PurpleCopyNumber> allCopyNumbers;
    private final AbstractFeatureReader<VariantContext, LineIterator> reader;
    private final StructuralVariantLegPloidyFactory<PurpleCopyNumber> ploidyFactory;

    public RecoverStructuralVariants(@NotNull final PurityAdjuster purityAdjuster, @NotNull final String recoveryVCF,
            @NotNull final List<PurpleCopyNumber> allCopyNumbers) {
        this.purityAdjuster = purityAdjuster;
        this.reader = getFeatureReader(recoveryVCF, new VCFCodec(), true);
        this.allCopyNumbers = Multimaps.fromRegions(allCopyNumbers);
        ploidyFactory = new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);
    }

    @NotNull
    public Collection<VariantContext> recoverVariants(@NotNull final List<StructuralVariant> currentVariants) throws IOException {
        final Map<String, VariantContext> result = Maps.newHashMap();

        final List<StructuralVariant> doubleEndedVariants =
                currentVariants.stream().filter(x -> x.type() != StructuralVariantType.SGL).collect(Collectors.toList());

        final List<StructuralVariantLegPloidy> doubleEndedPloidies = ploidyFactory.create(doubleEndedVariants, allCopyNumbers);
        final StructuralVariantLegCopyNumberChangeFactory changeFactory =
                new StructuralVariantLegCopyNumberChangeFactory(purityAdjuster, allCopyNumbers, currentVariants);

        recoverFromUnbalancedVariants(changeFactory, doubleEndedPloidies).forEach(x -> addToMap(result, x, "UNBALANCED_SV"));
        recoverFromUnexplainedSegments().forEach(x -> addToMap(result, x, "UNSUPPORTED_BREAKEND"));

        return result.values();
    }

    @NotNull
    private List<RecoveredVariant> recoverFromUnbalancedVariants(@NotNull final StructuralVariantLegCopyNumberChangeFactory changeFactory,
            @NotNull final List<StructuralVariantLegPloidy> svPloidies) throws IOException {
        final List<RecoveredVariant> result = Lists.newArrayList();

        for (StructuralVariantLegPloidy leg : svPloidies) {
            double copyNumber = leg.adjustedCopyNumber();
            double copyNumberChange = changeFactory.copyNumberChange(leg);
            double expectedCopyNumberChange = leg.averageImpliedPloidy();
            double unexplainedCopyNumberChange = Math.max(0, expectedCopyNumberChange - copyNumberChange);

            if (isUnbalanced(unexplainedCopyNumberChange, copyNumber)) {
                final List<PurpleCopyNumber> chromosomeCopyNumbers = allCopyNumbers.get(HumanChromosome.fromString(leg.chromosome()));
                int index = indexOf(leg.cnaPosition(), chromosomeCopyNumbers);
                if (index > 0) {
                    final PurpleCopyNumber prev = chromosomeCopyNumbers.get(index - 1);
                    final PurpleCopyNumber current = chromosomeCopyNumbers.get(index);

                    if (current.segmentStartSupport() != SegmentSupport.MULTIPLE && isSupportedByDepthWindowCounts(prev, current)) {
                        int expectedOrientation = -1 * leg.orientation();
                        final Optional<RecoveredVariant> recoveredVariant =
                                recoverSingleVariant(expectedOrientation, unexplainedCopyNumberChange, index, chromosomeCopyNumbers);
                        recoveredVariant.ifPresent(result::add);

                    }
                }
            }
        }

        return result;
    }

    @NotNull
    private List<RecoveredVariant> recoverFromUnexplainedSegments() throws IOException {

        final List<RecoveredVariant> result = Lists.newArrayList();

        for (Chromosome chromosome : allCopyNumbers.keySet()) {
            final List<PurpleCopyNumber> chromosomeCopyNumbers = allCopyNumbers.get(chromosome);

            for (int index = 1; index < chromosomeCopyNumbers.size() - 1; index++) {
                final PurpleCopyNumber current = chromosomeCopyNumbers.get(index);
                if (current.segmentStartSupport() == SegmentSupport.NONE) {
                    PurpleCopyNumber prev = chromosomeCopyNumbers.get(index - 1);
                    double unexplainedCopyNumberChange = Math.abs(prev.averageTumorCopyNumber() - current.averageTumorCopyNumber());

                    int expectedOrientation = Doubles.greaterThan(current.averageTumorCopyNumber(), prev.averageTumorCopyNumber()) ? -1 : 1;
                    recoverSingleVariant(expectedOrientation,
                            unexplainedCopyNumberChange,
                            index,
                            chromosomeCopyNumbers).ifPresent(result::add);
                }
            }
        }

        return result;
    }

    @NotNull
    private Optional<RecoveredVariant> recoverSingleVariant(int expectedOrientation, double unexplainedCopyNumberChange, int index,
            @NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        return recoverAllVariants(expectedOrientation, unexplainedCopyNumberChange, index, copyNumbers).stream().findFirst();
    }

    @NotNull
    private List<RecoveredVariant> recoverAllVariants(int expectedOrientation, double unexplainedCopyNumberChange, int index,
            @NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        assert (index > 1);

        PurpleCopyNumber prev = copyNumbers.get(index - 1);
        PurpleCopyNumber current = copyNumbers.get(index);

        final List<RecoveredVariant> result = Lists.newArrayList();
        long minPosition = Math.max(1, current.minStart() - 1000);
        long maxPosition = current.maxStart() + 1000;
        result.addAll(recover(expectedOrientation, unexplainedCopyNumberChange, minPosition, maxPosition, current, prev));
        result.sort(QUALITY_COMPARATOR.reversed());
        return result;
    }

    @NotNull
    private List<RecoveredVariant> recover(int expectedOrientation, double unexplainedCopyNumberChange, long min, long max,
            @NotNull final PurpleCopyNumber current, @NotNull final PurpleCopyNumber prev) throws IOException {
        final List<RecoveredVariant> result = Lists.newArrayList();

        final List<VariantContext> recovered = findVariants(current.chromosome(), min, max);
        for (VariantContext potentialVariant : recovered) {
            final String alt = potentialVariant.getAlternateAllele(0).getDisplayString();

            final String mateId = StructuralVariantFactory.mateId(potentialVariant);
            final String mateLocation = mateLocation(alt);
            final String mateChromosome = mateChromosome(mateLocation);
            final Long matePosition = matePosition(mateLocation);
            final int uncertainty = uncertainty(potentialVariant);

            final VariantContext mate = mateChromosome != null && matePosition != null && mateId != null ? findMate(mateId,
                    mateChromosome,
                    Math.max(1, matePosition - uncertainty),
                    matePosition + uncertainty) : null;

            final StructuralVariant sv = mate != null
                    ? StructuralVariantFactory.create(potentialVariant, mate)
                    : StructuralVariantFactory.createSingleBreakend(potentialVariant);

            final double ploidy = ploidyFactory.singleLegPloidy(sv.start(), prev.averageTumorCopyNumber(), current.averageTumorCopyNumber())
                    .averageImpliedPloidy();

            if (sv.start().orientation() == expectedOrientation && sufficientPloidy(ploidy, unexplainedCopyNumberChange)
                    && hasPotential(sv)) {
                result.add(ImmutableRecoveredVariant.builder()
                        .context(potentialVariant)
                        .mate(mate)
                        .variant(sv)
                        .copyNumber(current)
                        .prevCopyNumber(prev)
                        .build());
            }
        }

        return result;
    }

    private boolean sufficientPloidy(double ploidy, double unexplainedCopyNumberChange) {
        return Doubles.greaterOrEqual(ploidy, unexplainedCopyNumberChange * MIN_PLOIDY_AS_PERCENTAGE_OF_COPY_NUMBER_CHANGE)
                && Doubles.greaterOrEqual(ploidy, MIN_PLOIDY);
    }

    private int uncertainty(@NotNull final VariantContext context) {
        final int homlen = 2 * context.getAttributeAsInt("HOMLEN", 0);
        final int cipos = cipos(context);
        return Math.max(homlen, cipos);
    }

    private int cipos(@NotNull final VariantContext context) {
        int max = MIN_MATE_UNCERTAINTY;
        if (context.hasAttribute("IMPRECISE")) {

            final String cipos = context.getAttributeAsString("CIPOS", "-0,0");
            if (cipos.contains(",")) {
                for (String s : cipos.split(",")) {
                    try {
                        max = Math.max(max, 2 * Math.abs(Integer.valueOf(s)));
                    } catch (Exception ignored) {

                    }
                }
            }
        }
        return max;
    }

    private boolean hasPotential(@NotNull final StructuralVariant variant) {
        StructuralVariantLeg start = variant.start();
        StructuralVariantLeg end = variant.end();

        // This should never actually occur because we are searching within this area
        if (!isInRangeOfCopyNumberSegment(start, allCopyNumbers.get(HumanChromosome.fromString(start.chromosome())))) {
            return false;
        }

        if (end == null) {
            return variant.qualityScore() >= MIN_SINGLE_QUAL_SCORE;
        }

        if (variant.qualityScore() < MIN_MATE_QUAL_SCORE) {
            return false;
        }

        //        if (!isInRangeOfCopyNumberSegment(end, allCopyNumbers.get(HumanChromosome.fromString(end.chromosome())))) {
        //            return false;
        //        }

        long endPosition = end.position();
        StructuralVariantType type = variant.type();
        if (type == StructuralVariantType.DEL || type == StructuralVariantType.DUP || type == StructuralVariantType.INS) {
            assert (variant.end() != null);

            long length = Math.abs(endPosition - variant.start().position());
            return length >= MIN_LENGTH;
        }

        return true;
    }

    private boolean isInRangeOfCopyNumberSegment(@NotNull final StructuralVariantLeg leg,
            @NotNull final List<PurpleCopyNumber> copyNumbers) {
        final Predicate<PurpleCopyNumber> chrRange = copyNumber -> copyNumber.chromosome().equals(leg.chromosome());
        final Predicate<PurpleCopyNumber> posRange =
                copyNumber -> leg.cnaPosition() >= copyNumber.minStart() - 1000 && leg.cnaPosition() <= copyNumber.maxStart() + 1000;
        return copyNumbers.stream().anyMatch(chrRange.and(posRange));
    }

    @NotNull
    private List<VariantContext> findVariants(@NotNull final String chromosome, final long lowerBound, final long upperBound)
            throws IOException {

        return reader.query(chromosome, (int) lowerBound, (int) upperBound)
                .stream()
                .filter(RecoverStructuralVariants::isAppropriatelyFiltered)
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    static boolean isAppropriatelyFiltered(@NotNull VariantContext variantContext) {
        final Set<String> filters = variantContext.getFilters();
        return !filters.isEmpty() && !filters.contains(AF_FILTERED);
    }

    @NotNull
    private VariantContext findMate(@NotNull final String id, @NotNull final String chromosome, final long min, final long max)
            throws IOException {

        return reader.query(chromosome, (int) min, (int) max)
                .stream()
                .filter(x -> x.getID().equals(id))
                .findFirst()
                .orElseThrow(() -> new IOException("Unable to find mateId " + id + " between " + min + " and " + max));
    }

    @Nullable
    static String mateLocation(@NotNull final String alt) {
        final String bracket;
        if (alt.contains("[")) {
            bracket = "\\[";
        } else if (alt.contains("]")) {
            bracket = "]";
        } else {
            return null;
        }

        String[] results = alt.split(bracket);
        for (String result : results) {
            if (result.contains(":")) {
                return result;
            }
        }
        return null;
    }

    @Nullable
    private static String mateChromosome(@Nullable String mate) {
        return mate == null || !mate.contains(":") ? null : mate.split(":")[0];
    }

    @Nullable
    private static Long matePosition(@Nullable String mate) {
        return mate == null || !mate.contains(":") ? null : Long.valueOf(mate.split(":")[1]);
    }

    @VisibleForTesting
    private static <T extends GenomeRegion> int indexOf(long cnaPosition, @NotNull final List<T> regions) {
        assert (!regions.isEmpty());
        for (int i = 0; i < regions.size(); i++) {
            if (regions.get(i).start() == cnaPosition) {
                return i;
            }
        }

        return -1;
    }

    @VisibleForTesting
    @NotNull
    static <T extends GenomeRegion> T closest(long position, @NotNull final List<T> regions) {
        assert (!regions.isEmpty());

        long minDistance = position - 1;
        for (int i = 1; i < regions.size(); i++) {
            long distanceFromRegion = position - regions.get(i).start();
            if (distanceFromRegion < 0) {
                return Math.abs(distanceFromRegion) < minDistance ? regions.get(i) : regions.get(i - 1);
            }

            minDistance = distanceFromRegion;
        }

        return regions.get(regions.size() - 1);
    }

    private static void addToMap(@NotNull final Map<String, VariantContext> map, @NotNull final RecoveredVariant variant,
            @NotNull final String method) {
        recover(variant, method).forEach(x -> map.put(x.getID(), x));
    }

    @NotNull
    private static List<VariantContext> recover(@NotNull final RecoveredVariant variant, @NotNull final String partialMethod) {
        final List<VariantContext> result = Lists.newArrayList();
        final Set<String> recoveryFilterSet = filterSet(variant.context());
        final VariantContext mate = variant.mate();

        final String recoveryMethod;
        final String recoveryFilter;
        if (mate != null) {
            final GenomePosition matePosition = GenomePositions.create(mate.getContig(), mate.getStart());
            final GenomePosition contextPosition = GenomePositions.create(variant.context().getContig(), variant.context().getStart());
            final boolean contextIsStart = contextPosition.compareTo(matePosition) <= 0;
            recoveryFilterSet.addAll(filterSet(mate));

            recoveryMethod = partialMethod + (contextIsStart ? "_START" : "_END");
            recoveryFilter = filterString(recoveryFilterSet);
            result.add(recover(mate, recoveryMethod, recoveryFilter));

        } else {
            recoveryMethod = partialMethod + "_START";
            recoveryFilter = filterString(recoveryFilterSet);
        }

        result.add(recover(variant.context(), recoveryMethod, recoveryFilter));
        return result;
    }

    @NotNull
    private static String filterString(@NotNull final Set<String> filters) {
        if (filters.isEmpty()) {
            return "PASS";
        }

        final StringJoiner recoveryFilterJoiner = new StringJoiner(",");
        filters.stream().sorted().forEach(recoveryFilterJoiner::add);

        return recoveryFilterJoiner.toString();
    }

    @NotNull
    private static Set<String> filterSet(@NotNull VariantContext variantContext) {
        return variantContext.isNotFiltered() ? Sets.newHashSet("PASS") : Sets.newHashSet(variantContext.getFilters());
    }

    @NotNull
    private static VariantContext recover(@NotNull final VariantContext context, @NotNull final String recoveryMethod,
            @NotNull final String recoveryFilter) {
        return new VariantContextBuilder(context).unfiltered()
                .attribute(StructuralVariantFactory.RECOVERED, true)
                .attribute(RECOVERY_METHOD, recoveryMethod)
                .attribute(RECOVERY_FILTER, recoveryFilter)
                .make();
    }

    private static boolean isUnbalanced(double unexplainedCopyNumberChange, double copyNumber) {
        return Doubles.greaterOrEqual(unexplainedCopyNumberChange,
                UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE_AS_PERCENT_OF_COPY_NUMBER * copyNumber) && Doubles.greaterOrEqual(
                unexplainedCopyNumberChange,
                UNBALANCED_MIN_UNEXPLAINED_COPY_NUMBER_CHANGE);
    }

    private static boolean isSupportedByDepthWindowCounts(@NotNull final PurpleCopyNumber prev, @Nullable final PurpleCopyNumber next) {
        return prev.depthWindowCount() >= UNBALANCED_MIN_DEPTH_WINDOW_COUNT && (next == null
                || next.depthWindowCount() >= UNBALANCED_MIN_DEPTH_WINDOW_COUNT);
    }

    @Override
    public void close() throws IOException {
        reader.close();
    }
}

package com.hartwig.hmftools.common.variant.recovery;

import static java.util.Comparator.comparingDouble;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.Closeable;
import java.io.IOException;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
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
import htsjdk.variant.vcf.VCFCodec;

public class RecoverStructuralVariants implements Closeable {

    private static final double MIN_LENGTH = 1000;
    private static final double MIN_MATE_QUAL_SCORE = 350;
    private static final double MIN_SINGLE_QUAL_SCORE = 1000;
    private static final double UNBALANCED_MIN_PLOIDY = 0.5;
    private static final double UNBALANCED_MAX_COPY_NUMBER_CHANGE = 0.25;
    private static final double UNBALANCED_MIN_DEPTH_WINDOW_COUNT = 5;
    private static final int MIN_MATE_UNCERTAINTY = 150;
    private static final String AF_FILTERED = "af";
    private static final String PON_FILTERED = "PON";

    private static final Comparator<RecoveredVariant> QUALITY_COMPARATOR = comparingDouble(x -> x.context().getPhredScaledQual());

    private final ListMultimap<String, PurpleCopyNumber> allCopyNumbers;
    private final AbstractFeatureReader<VariantContext, LineIterator> reader;
    private final StructuralVariantLegPloidyFactory<PurpleCopyNumber> ploidyFactory;

    public RecoverStructuralVariants(@NotNull final PurityAdjuster purityAdjuster, @NotNull final String recoveryVCF,
            @NotNull final List<PurpleCopyNumber> allCopyNumbers) {
        this.reader = getFeatureReader(recoveryVCF, new VCFCodec(), true);
        this.allCopyNumbers = ArrayListMultimap.create();
        for (PurpleCopyNumber copyNumber : allCopyNumbers) {
            this.allCopyNumbers.put(copyNumber.chromosome(), copyNumber);
        }
        ploidyFactory = new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);
    }

    @NotNull
    public Collection<VariantContext> recoverVariants(@NotNull final List<StructuralVariant> currentVariants) throws IOException {
        final Map<String, VariantContext> result = Maps.newHashMap();

        final List<StructuralVariant> doubleEndedVariants =
                currentVariants.stream().filter(x -> x.type() != StructuralVariantType.SGL).collect(Collectors.toList());

        final List<StructuralVariantLegPloidy> doubleEndedPloidies = ploidyFactory.create(doubleEndedVariants, allCopyNumbers);

        recoverFromUnbalancedVariants(doubleEndedPloidies).forEach(x -> addToMap(result, x));
        recoverFromUnexplainedSegments().forEach(x -> addToMap(result, x));

        return result.values();
    }

    @NotNull
    private List<RecoveredVariant> recoverFromUnbalancedVariants(@NotNull final List<StructuralVariantLegPloidy> svPloidies)
            throws IOException {
        final List<RecoveredVariant> result = Lists.newArrayList();

        for (StructuralVariantLegPloidy svPloidy : svPloidies) {
            if (isUnbalanced(svPloidy)) {
                final List<PurpleCopyNumber> chromosomeCopyNumbers = allCopyNumbers.get(svPloidy.chromosome());
                int index = indexOf(svPloidy.position(), chromosomeCopyNumbers);
                if (index > 1) {
                    final PurpleCopyNumber prev = chromosomeCopyNumbers.get(index - 1);
                    final PurpleCopyNumber next = index < chromosomeCopyNumbers.size() - 2 ? chromosomeCopyNumbers.get(index + 1) : null;
                    if (isSupportedByDepthWindowCounts(prev, next)) {
                        int expectedOrientation = -1 * svPloidy.orientation();
                        recoverSingleVariant(expectedOrientation, index, chromosomeCopyNumbers).ifPresent(result::add);
                    }
                }
            }
        }

        return result;
    }

    @NotNull
    private List<RecoveredVariant> recoverFromUnexplainedSegments() throws IOException {

        final List<RecoveredVariant> result = Lists.newArrayList();

        for (String chromosome : allCopyNumbers.keySet()) {
            final List<PurpleCopyNumber> chromosomeCopyNumbers = allCopyNumbers.get(chromosome);

            for (int index = 1; index < chromosomeCopyNumbers.size() - 1; index++) {
                final PurpleCopyNumber current = chromosomeCopyNumbers.get(index);
                if (current.segmentStartSupport() == SegmentSupport.NONE) {
                    PurpleCopyNumber prev = chromosomeCopyNumbers.get(index - 1);
                    int expectedOrientation = Doubles.greaterThan(current.averageTumorCopyNumber(), prev.averageTumorCopyNumber()) ? -1 : 1;
                    recoverSingleVariant(expectedOrientation, index, chromosomeCopyNumbers).ifPresent(result::add);
                }
            }
        }

        return result;
    }

    @NotNull
    private Optional<RecoveredVariant> recoverSingleVariant(int expectedOrientation, int index,
            @NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        return recoverAllVariants(expectedOrientation, index, copyNumbers).stream().findFirst();
    }

    @NotNull
    private List<RecoveredVariant> recoverAllVariants(int expectedOrientation, int index, @NotNull final List<PurpleCopyNumber> copyNumbers)
            throws IOException {
        assert (index > 1);

        PurpleCopyNumber prev = copyNumbers.get(index - 1);
        PurpleCopyNumber current = copyNumbers.get(index);

        final List<RecoveredVariant> result = Lists.newArrayList();
        long minPosition = Math.max(1, current.minStart() - 1000);
        long maxPosition = current.maxStart() + 1000;
        result.addAll(recover(expectedOrientation, minPosition, maxPosition, current, prev));
        result.sort(QUALITY_COMPARATOR.reversed());
        return result;
    }

    @NotNull
    private List<RecoveredVariant> recover(int expectedOrientation, long min, long max, @NotNull final PurpleCopyNumber current,
            @NotNull final PurpleCopyNumber prev) throws IOException {
        List<RecoveredVariant> result = Lists.newArrayList();

        List<VariantContext> recovered = findVariants(current.chromosome(), min, max);
        for (VariantContext potentialVariant : recovered) {
            final String alt = potentialVariant.getAlternateAllele(0).getDisplayString();
            int orientation = orientation(alt);
            if (orientation != expectedOrientation) {
                continue;
            }

            final String mateId = StructuralVariantFactory.mateId(potentialVariant);
            final String mateLocation = mateLocation(alt);
            final String mateChromosome = mateChromosome(mateLocation);
            final Long matePosition = matePosition(mateLocation);
            final int uncertainty = uncertainty(potentialVariant);

            final VariantContext mate = mateChromosome != null && matePosition != null && mateId != null ? findMate(mateId,
                    mateChromosome,
                    Math.max(1, matePosition - uncertainty),
                    matePosition + uncertainty) : null;

            final PurpleCopyNumber mateCopyNumber = mate == null || !allCopyNumbers.containsKey(mate.getContig())
                    ? null
                    : closest(mate.getStart(), allCopyNumbers.get(mate.getContig()));

            final StructuralVariant sv = mate != null
                    ? StructuralVariantFactory.create(potentialVariant, mate)
                    : StructuralVariantFactory.createSingleBreakend(potentialVariant);

            if (hasPotential(min, max, sv)) {
                result.add(ImmutableRecoveredVariant.builder()
                        .context(potentialVariant)
                        .mate(mate)
                        .mateCopyNumber(mateCopyNumber)
                        .variant(sv)
                        .copyNumber(current)
                        .prevCopyNumber(prev)
                        .build());
            }
        }

        return result;
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

    private boolean hasPotential(long min, long max, @NotNull StructuralVariant variant) {
        StructuralVariantLeg end = variant.end();
        if (end == null) {
            return variant.qualityScore() >= MIN_SINGLE_QUAL_SCORE;
        }

        if (variant.qualityScore() < MIN_MATE_QUAL_SCORE) {
            return false;
        }

        long endPosition = end.position();
        StructuralVariantType type = variant.type();
        if (type == StructuralVariantType.DEL || type == StructuralVariantType.DUP || type == StructuralVariantType.INS) {
            assert (variant.end() != null);

            long length = Math.abs(endPosition - variant.start().position());
            if (length < MIN_LENGTH) {
                return false;
            }

            long variantStart = Math.min(variant.start().position(), endPosition);
            long variantEnd = Math.max(variant.start().position(), endPosition);
            return variantStart <= min || variantEnd >= max;
        }

        return true;
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
        return !filters.isEmpty() && !filters.contains(PON_FILTERED) && !filters.contains(AF_FILTERED);
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
    static int orientation(@NotNull final String alt) {
        if (alt.charAt(0) == '.') {
            return -1;
        }

        if (alt.charAt(alt.length() - 1) == '.') {
            return 1;
        }

        if (alt.contains("[")) {
            return 1;
        }

        if (alt.contains("]")) {
            return -1;
        }

        return 0;
    }

    @VisibleForTesting
    private static <T extends GenomeRegion> int indexOf(long start, @NotNull final List<T> regions) {
        assert (!regions.isEmpty());
        for (int i = 0; i < regions.size(); i++) {
            if (regions.get(i).start() == start) {
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

    private static void addToMap(@NotNull Map<String, VariantContext> map, @NotNull RecoveredVariant variant) {
        map.put(variant.context().getID(), variant.context());
        VariantContext mate = variant.mate();
        if (mate != null) {
            map.put(mate.getID(), mate);
        }
    }

    private static boolean isUnbalanced(@NotNull final StructuralVariantLegPloidy svPloidy) {
        return Doubles.greaterThan(svPloidy.averageImpliedPloidy(), UNBALANCED_MIN_PLOIDY) && Doubles.lessThan(absAdjustedCopyNumberChange(
                svPloidy), UNBALANCED_MAX_COPY_NUMBER_CHANGE);
    }

    private static boolean isSupportedByDepthWindowCounts(@NotNull final PurpleCopyNumber prev, @Nullable final PurpleCopyNumber next) {
        return prev.depthWindowCount() >= UNBALANCED_MIN_DEPTH_WINDOW_COUNT && (next == null || next.depthWindowCount()
                >= UNBALANCED_MIN_DEPTH_WINDOW_COUNT);
    }

    private static double absAdjustedCopyNumberChange(@NotNull final StructuralVariantLegPloidy ploidy) {
        double leftCopyNumber = ploidy.leftCopyNumber().orElse(0D);
        double rightCopyNumber = ploidy.rightCopyNumber().orElse(0D);

        return ploidy.orientation() == 1 ? leftCopyNumber - rightCopyNumber : rightCopyNumber - leftCopyNumber;
    }

    @Override
    public void close() throws IOException {
        reader.close();
    }
}

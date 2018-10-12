package com.hartwig.hmftools.common.variant.recovery;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class StructuralVariantRecovery {

    private static final double MIN_SINGLE_QUAL_SCORE = 250; //TODO, FIX THIS AGAIN
    private static final double MIN_MATE_QUAL_SCORE = 250;
    private static final double MIN_LENGTH = 1000;

    private final AbstractFeatureReader<VariantContext, LineIterator> reader;
    private final ListMultimap<String, PurpleCopyNumber> allCopyNumbers;

    public StructuralVariantRecovery(@NotNull final String vcfFile, @NotNull final ListMultimap<String, PurpleCopyNumber> allCopyNumbers) {
        this.reader = getFeatureReader(vcfFile, new VCFCodec(), true);
        this.allCopyNumbers = allCopyNumbers;
    }

    public List<RecoveredVariant> recoverUnbalancedVariants(List<StructuralVariantLegPloidy> svPloidies) throws IOException {
        final List<RecoveredVariant> result = Lists.newArrayList();

        for (StructuralVariantLegPloidy svPloidy : svPloidies) {
            if (Doubles.greaterThan(svPloidy.averageImpliedPloidy(), 0.5) && absAdjustedCopyNumberChange(svPloidy) < 0.25) {
                final List<PurpleCopyNumber> chromosomeCopyNumbers = allCopyNumbers.get(svPloidy.chromosome());
                int index = indexOf(svPloidy.position(), chromosomeCopyNumbers);
                if (index > 1) {
                    result.addAll(recoverVariants(index, chromosomeCopyNumbers));
                }
            }
        }

        return result;
    }

    private static double absAdjustedCopyNumberChange(@NotNull final StructuralVariantLegPloidy ploidy) {
        double leftCopyNumber = ploidy.leftCopyNumber().orElse(0D);
        double rightCopyNumber = ploidy.rightCopyNumber().orElse(0D);

        return ploidy.orientation() == 1 ? leftCopyNumber - rightCopyNumber : rightCopyNumber - leftCopyNumber;
    }

    public List<RecoveredVariant> recoverVariants() throws IOException {

        final List<RecoveredVariant> result = Lists.newArrayList();

        for (String chromosome : allCopyNumbers.keySet()) {
            final List<PurpleCopyNumber> chromosomeCopyNumbers = allCopyNumbers.get(chromosome);

            for (int i = 1; i < chromosomeCopyNumbers.size() - 1; i++) {
                final PurpleCopyNumber current = chromosomeCopyNumbers.get(i);
                if (current.segmentStartSupport() == SegmentSupport.NONE) {
                    result.addAll(recoverVariants(i, chromosomeCopyNumbers));
                }
            }
        }

        return result;
    }

    @NotNull
    private List<RecoveredVariant> recoverVariants(int index, @NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {
        assert (index > 1);

        PurpleCopyNumber prev = copyNumbers.get(index - 1);
        PurpleCopyNumber current = copyNumbers.get(index);
        PurpleCopyNumber next = index < copyNumbers.size() - 1 ? copyNumbers.get(index + 1) : null;

        final List<RecoveredVariant> result = Lists.newArrayList();
        long minPosition = Math.max(1, current.minStart() - 1000);
        long maxPosition = current.maxStart() + 1000;
        result.addAll(recover(minPosition, maxPosition, current, prev, next, allCopyNumbers));
        return result;
    }

    @NotNull
    public List<RecoveredVariant> recover(long min, long max, @NotNull final PurpleCopyNumber current, @NotNull final PurpleCopyNumber prev,
            @Nullable final PurpleCopyNumber next, @NotNull final ListMultimap<String, PurpleCopyNumber> allCopyNumbers)
            throws IOException {

        List<RecoveredVariant> result = Lists.newArrayList();

        ImmutableRecoveredVariant.Builder builder = ImmutableRecoveredVariant.builder()
                .from(current)
                .minStart(current.minStart())
                .maxStart(current.maxStart())
                .copyNumber(current.averageTumorCopyNumber())
                .support(current.segmentStartSupport())
                .next(current.segmentEndSupport())
                .previous(prev.segmentStartSupport())
                .baf(current.averageActualBAF())
                .depthWindowCount(current.depthWindowCount())
                .gcContent(current.gcContent())
                .prevGCContent(prev.gcContent())
                .nextGCContent(next == null ? -1 : next.gcContent())
                .prevLength(prev.end() - prev.start() + 1)
                .prevCopyNumber(prev.averageTumorCopyNumber())
                .prevBaf(prev.averageActualBAF())
                .prevDepthWindowCount(prev.depthWindowCount());

        int expectedOrientation = Doubles.greaterThan(current.averageTumorCopyNumber(), prev.averageTumorCopyNumber()) ? -1 : 1;

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

            final StructuralVariant sv = mate != null
                    ? StructuralVariantFactory.create(potentialVariant, mate)
                    : StructuralVariantFactory.createSingleBreakend(potentialVariant);

            if (hasPotential(min, max, sv)) {
                final StructuralVariantLeg end = sv.end();
                final PurpleCopyNumber mateCopyNumber = end == null || !allCopyNumbers.containsKey(end.chromosome()) ? null : closest(end.position(), allCopyNumbers.get(end.chromosome()));

                builder.alt(alt)
                        .qual(sv.qualityScore())
                        .variant(sv.start().chromosome() + ":" + sv.start().position())
                        .orientation(orientation)
                        .tumourReferenceFragmentCount(sv.start().tumourReferenceFragmentCount())
                        .tumourVariantFragmentCount(sv.start().tumourVariantFragmentCount())
                        .mate(end == null ? null : end.chromosome() + ":" + end.position())
                        .mateOrientation(end == null ? null : (int) end.orientation())
                        .mateMinStart(mateCopyNumber == null ? null : mateCopyNumber.minStart())
                        .mateMaxStart(mateCopyNumber == null ? null : mateCopyNumber.maxStart())
                        .mateSupport(mateCopyNumber == null ? null : mateCopyNumber.segmentStartSupport())
                        .mateTumourReferenceFragmentCount(end == null ? null : end.tumourReferenceFragmentCount())
                        .mateTumourVariantFragmentCount(end == null ? null : end.tumourVariantFragmentCount())
                        .filter(sv.filter());

                result.add(builder.build());
            }
        }

        if (result.isEmpty()) {
            result.add(builder.build());
        }

        return result;
    }

    private int uncertainty(@NotNull final VariantContext context) {
        final int homlen = 2 * context.getAttributeAsInt("HOMLEN", 0);
        final int cipos = cipos(context);
        return Math.max(homlen, cipos);
    }

    private int cipos(@NotNull final VariantContext context) {
        int max = 150;
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
        final List<VariantContext> result = Lists.newArrayList();

        try (CloseableTribbleIterator<VariantContext> iterator = reader.query(chromosome, (int) lowerBound, (int) upperBound)) {
            for (VariantContext variant : iterator) {
                if (variant.isFiltered()) {
                    result.add(variant);
                }
            }
        }

        return result;
    }

    @NotNull
    private VariantContext findMate(@NotNull final String id, @NotNull final String chromosome, final long min, final long max)
            throws IOException {

        try (CloseableTribbleIterator<VariantContext> iterator = reader.query(chromosome, (int) min, (int) max)) {
            for (VariantContext variant : iterator) {
                if (variant.getID().equals(id)) {
                    return variant;
                }
            }
        }

        throw new IOException("Unable to find mateId " + id + " between " + min + " and " + max);
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
    static String mateChromosome(@Nullable String mate) {
        return mate == null || !mate.contains(":") ? null : mate.split(":")[0];
    }

    @Nullable
    static Long matePosition(@Nullable String mate) {
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
    static <T extends GenomeRegion> int indexOf(long start, @NotNull final List<T> regions) {
        assert (!regions.isEmpty());
        for (int i = 0; i < regions.size(); i++) {
            if (regions.get(i).start() == start) {
                return i;
            }
        }

        throw new UnsupportedOperationException();

    }

    @VisibleForTesting
    @NotNull
    static <T extends GenomeRegion> T closest(long position, @NotNull final List<T> regions) {
        assert(!regions.isEmpty());


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

}

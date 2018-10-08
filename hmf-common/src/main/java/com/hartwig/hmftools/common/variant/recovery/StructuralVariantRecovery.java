package com.hartwig.hmftools.common.variant.recovery;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
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

    private static final double MIN_SINGLE_QUAL_SCORE = 1000;
    private static final double MIN_MATE_QUAL_SCORE = 250;
    private static final double MIN_LENGTH = 1000;

    private final AbstractFeatureReader<VariantContext, LineIterator> reader;

    public StructuralVariantRecovery(@NotNull final String vcfFile) {
        this.reader = getFeatureReader(vcfFile, new VCFCodec(), true);
    }

    public void doStuff2(List<StructuralVariantLegPloidy> svPloidies) throws IOException {
        for (StructuralVariantLegPloidy svPloidy : svPloidies) {
            if (Doubles.greaterThan(svPloidy.averageImpliedPloidy(), 0.5) && adjustedCopyNumberChange(svPloidy) > 0.15) {

                System.out.println();

                System.out.println("MADE IT!!");
            }
        }
    }

    @VisibleForTesting
    static double adjustedCopyNumberChange(@NotNull final StructuralVariantLegPloidy ploidy) {
        double leftCopyNumber = ploidy.leftCopyNumber().orElse(0D);
        double rightCopyNumber = ploidy.rightCopyNumber().orElse(0D);

        return ploidy.orientation() == 1 ? leftCopyNumber - rightCopyNumber : rightCopyNumber - leftCopyNumber;
    }

    public List<RecoveredVariant> doStuff(@NotNull final List<PurpleCopyNumber> copyNumbers) throws IOException {

        final List<RecoveredVariant> result = Lists.newArrayList();

        for (int i = 1; i < copyNumbers.size() - 1; i++) {
            PurpleCopyNumber prev = copyNumbers.get(i - 1);
            PurpleCopyNumber current = copyNumbers.get(i);
            PurpleCopyNumber next = copyNumbers.get(i + 1);

            if (current.segmentStartSupport() == SegmentSupport.NONE) {
                long minPosition = current.minStart() - 1000;
                long maxPosition = current.maxStart() + 1000;
                result.addAll(recover(minPosition, maxPosition, current, prev, next));
            }
        }

        return result;
    }

    @NotNull
    public List<RecoveredVariant> recover(long min, long max, @NotNull final PurpleCopyNumber current, @NotNull final PurpleCopyNumber prev,
            @NotNull final PurpleCopyNumber next) throws IOException {
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
                .nextGCContent(next.gcContent())
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
            final Optional<VariantContext> mate;
            if (mateChromosome != null && matePosition != null && mateId != null) {
                mate = findMate(mateId, mateChromosome, matePosition);
            } else {
                mate = Optional.empty();
            }

            if (mate.isPresent()) {
                StructuralVariant sv = StructuralVariantFactory.create(potentialVariant, mate.get());
                if (hasPotential(min, max, sv)) {
                    builder.alt(alt)
                            .qual(potentialVariant.getPhredScaledQual())
                            .variant(potentialVariant.getContig() + ":" + potentialVariant.getStart())
                            .orientation(orientation)
                            .mate(sv.end().chromosome() + ":" + sv.end().position())
                            .mateOrientation((int) sv.end().orientation())
                            .filter(filter(potentialVariant.getFilters()));
                }
            } else {
                if (Doubles.greaterOrEqual(potentialVariant.getPhredScaledQual(), MIN_SINGLE_QUAL_SCORE)) {
                    builder.alt(alt)
                            .qual(potentialVariant.getPhredScaledQual())
                            .variant(potentialVariant.getContig() + ":" + potentialVariant.getStart())
                            .orientation(orientation)
                            .mate(null)
                            .mateOrientation(null)
                            .filter(filter(potentialVariant.getFilters()));
                    result.add(builder.build());
                }
            }

        }

        if (result.isEmpty()) {
            result.add(builder.build());
        }

        return result;
    }

    private boolean hasPotential(long min, long max, @NotNull StructuralVariant variant) {
        StructuralVariantLeg end = variant.end();
        assert (end != null);

        if (variant.qualityScore() < MIN_MATE_QUAL_SCORE) {
            return false;
        }

        long endPosition = end.position();
        StructuralVariantType type = variant.type();
        if (type == StructuralVariantType.DEL || type == StructuralVariantType.DUP) {
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
                result.add(variant);
            }
        }

        return result;
    }

    @NotNull
    private Optional<VariantContext> findMate(@NotNull final String id, @NotNull final String chromosome, final long position)
            throws IOException {

        try (CloseableTribbleIterator<VariantContext> iterator = reader.query(chromosome, (int) position, (int) position)) {
            for (VariantContext variant : iterator) {
                if (variant.getID().equals(id)) {
                    return Optional.of(variant);
                }
            }
        }

        return Optional.empty();
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

    @NotNull
    private static String filter(@NotNull final Set<String> filters) {
        if (filters.isEmpty()) {
            return "PASS";
        }

        StringJoiner joiner = new StringJoiner(",");
        filters.forEach(joiner::add);
        return joiner.toString();
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

}

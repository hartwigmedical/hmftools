package com.hartwig.hmftools.common.variant.recovery;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class StructuralVariantRecovery {

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
                long minPosition = Math.max(prev.end() + 1, current.minStart() - 1000);
                long maxPosition = Math.min(next.start() - 1, current.maxStart() + 1000);
                result.addAll(recover(minPosition, maxPosition, current, prev, next));

            }
        }

        return result;
    }

    @NotNull
    public List<RecoveredVariant> recover(long min, long max, @NotNull final PurpleCopyNumber current, @NotNull final PurpleCopyNumber prev, @NotNull final PurpleCopyNumber next)
            throws IOException {
        List<RecoveredVariant> result = Lists.newArrayList();

        ImmutableRecoveredVariant.Builder builder = ImmutableRecoveredVariant.builder()
                .from(current)
                .copyNumber(current.averageTumorCopyNumber())
                .support(current.segmentStartSupport())
                .next(current.segmentEndSupport())
                .previous(prev.segmentStartSupport())
                .baf(current.averageActualBAF())
                .depthWindowCount(current.depthWindowCount())
                .gcContent(current.gcContent())
                .prevGCContent(prev.gcContent())
                .nextGCContent(next.gcContent())
                .prevLength(prev.end() - prev.start()  + 1)
                .prevCopyNumber(prev.averageTumorCopyNumber())
                .prevBaf(prev.averageActualBAF())
                .prevDepthWindowCount(prev.depthWindowCount());

        List<VariantContext> recovered = findVariants(current.chromosome(), min, max);
        for (VariantContext potentialVariant : recovered) {

            final String alt = potentialVariant.getAlternateAllele(0).getDisplayString();
            builder.alt(alt)
                    .qual(potentialVariant.getPhredScaledQual())
                    .variant(potentialVariant.getContig() + ":" + potentialVariant.getStart())
                    .orientation(orientation(alt))
                    .mate(mate(alt))
                    .filter(filter(potentialVariant.getFilters()));

            result.add(builder.build());

        }

        if (result.isEmpty()) {
            result.add(builder.build());
        }

        return result;
    }

    @NotNull
    public List<VariantContext> findVariants(@NotNull final String chromosome, final long lowerBound, final long upperBound)
            throws IOException {
        final List<VariantContext> result = Lists.newArrayList();

        try (CloseableTribbleIterator<VariantContext> iterator = reader.query(chromosome, (int) lowerBound, (int) upperBound)) {
            for (VariantContext variant : iterator) {
                result.add(variant);
            }
        }

        return result;
    }

    @Nullable
    private static String mate(@NotNull final String alt) {
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

    @NotNull
    private static String filter(@NotNull final Set<String> filters) {
        if (filters.isEmpty()) {
            return "PASS";
        }

        StringJoiner joiner = new StringJoiner(",");
        filters.forEach(joiner::add);
        return joiner.toString();
    }

    private static int orientation(@NotNull final String alt) {
        if (alt.charAt(0) == '.') {
            return -1;
        }

        if (alt.charAt(alt.length() - 1) == '.') {
            return 1;
        }

        if (alt.contains("[")) {
            return 1;
        }

        if (alt.contains("\\]")) {
            return -1;
        }

        return 0;
    }

}

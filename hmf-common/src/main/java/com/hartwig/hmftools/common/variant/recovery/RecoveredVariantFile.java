package com.hartwig.hmftools.common.variant.recovery;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantLeg;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class RecoveredVariantFile {

    private static final String HEADER_PREFIX = "#";
    private static final String DELIMITER = "\t";
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00");

    public static void write(@NotNull final String filePath, @NotNull Collection<RecoveredContext> variants) throws IOException {
        final Collection<String> lines = Lists.newArrayList();
        lines.add(header());
        variants.stream().map(RecoveredVariantFile::toString).forEach(lines::add);
        Files.write(new File(filePath).toPath(), lines);
    }

    @NotNull
    public static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "").add("chromosome")
                .add("start")
                .add("minStart")
                .add("maxStart")
                .add("end")
                .add("bases")
                .add("baf")
                .add("copyNumber")
                .add("depthWindowCount")
                .add("gcContent")
                .add("support")
                .add("tumourVariantFragmentCount")
                .add("tumourReferenceFragmentCount")
                .add("prevBases")
                .add("prevBaf")
                .add("prevCopyNumber")
                .add("prevDepthWindowCount")
                .add("prevGCContent")
                .add("previousSupport")
                .add("nextGCContent")
                .add("nextSupport")
                .add("variant")
                .add("orientation")
                .add("qual")
                .add("filter")
                .add("mate")
                .add("mateOrientation")
                .add("mateMinStart")
                .add("mateMaxStart")
                .add("mateSupport")
                .add("mateTumourVariantFragmentCount")
                .add("mateTumourReferenceFragmentCount")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final RecoveredContext recoveredVariant) {

        StructuralVariant variant = recoveredVariant.variant();
        StructuralVariantLeg start = variant.start();
        @Nullable
        StructuralVariantLeg end = variant.end();

        PurpleCopyNumber copyNumber = recoveredVariant.copyNumber();
        PurpleCopyNumber preCopyNumber = recoveredVariant.prevCopyNumber();
        PurpleCopyNumber mateCopyNumber = recoveredVariant.mateCopyNumber();

        return new StringJoiner(DELIMITER).add(String.valueOf(copyNumber.chromosome()))
                .add(String.valueOf(copyNumber.start()))
                .add(String.valueOf(copyNumber.minStart()))
                .add(String.valueOf(copyNumber.maxStart()))
                .add(String.valueOf(copyNumber.end()))
                .add(String.valueOf(copyNumber.bases()))
                .add(FORMAT.format(copyNumber.averageActualBAF()))
                .add(FORMAT.format(copyNumber.averageTumorCopyNumber()))
                .add(String.valueOf(copyNumber.depthWindowCount()))
                .add(String.valueOf(copyNumber.gcContent()))
                .add(String.valueOf(copyNumber.segmentStartSupport()))
                .add(String.valueOf(variant.start().tumourVariantFragmentCount()))
                .add(String.valueOf(variant.start().tumourReferenceFragmentCount()))

                .add(String.valueOf(preCopyNumber.bases()))
                .add(FORMAT.format(preCopyNumber.averageActualBAF()))
                .add(FORMAT.format(preCopyNumber.averageTumorCopyNumber()))
                .add(String.valueOf(preCopyNumber.depthWindowCount()))
                .add(String.valueOf(preCopyNumber.gcContent()))
                .add(String.valueOf(preCopyNumber.segmentStartSupport()))

                .add("-1")
                .add(String.valueOf(copyNumber.segmentEndSupport()))
                .add(String.valueOf(start.chromosome() + ":" + start.position()))
                .add(String.valueOf(start.orientation()))
                .add(String.valueOf(variant.qualityScore()))
                .add(String.valueOf(variant.filter()))
                .add(String.valueOf(recoveredVariant.mate()))
                .add(String.valueOf(end == null ? null : end.orientation()))
                .add(String.valueOf(mateCopyNumber == null ? null : mateCopyNumber.minStart()))
                .add(String.valueOf(mateCopyNumber == null ? null : mateCopyNumber.maxStart()))
                .add(String.valueOf(mateCopyNumber == null ? null : mateCopyNumber.segmentStartSupport()))
                .add(String.valueOf(end == null ? null : end.tumourVariantFragmentCount()))
                .add(String.valueOf(end == null ? null : end.tumourReferenceFragmentCount()))
                .toString();
    }

}

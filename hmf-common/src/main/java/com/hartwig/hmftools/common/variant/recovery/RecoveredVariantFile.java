package com.hartwig.hmftools.common.variant.recovery;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class RecoveredVariantFile {

    private static final String HEADER_PREFIX = "#";
    private static final String DELIMITER = "\t";
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00");

    public static void write(@NotNull final String filePath, @NotNull Collection<RecoveredVariant> variants) throws IOException {
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
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final RecoveredVariant recoveredVariant) {
        return new StringJoiner(DELIMITER).add(String.valueOf(recoveredVariant.chromosome()))
                .add(String.valueOf(recoveredVariant.start()))
                .add(String.valueOf(recoveredVariant.minStart()))
                .add(String.valueOf(recoveredVariant.maxStart()))
                .add(String.valueOf(recoveredVariant.end()))
                .add(String.valueOf(recoveredVariant.bases()))
                .add(FORMAT.format(recoveredVariant.baf()))
                .add(FORMAT.format(recoveredVariant.copyNumber()))
                .add(String.valueOf(recoveredVariant.depthWindowCount()))
                .add(String.valueOf(recoveredVariant.gcContent()))

                .add(String.valueOf(recoveredVariant.prevLength()))
                .add(FORMAT.format(recoveredVariant.prevBaf()))
                .add(FORMAT.format(recoveredVariant.prevCopyNumber()))
                .add(String.valueOf(recoveredVariant.prevDepthWindowCount()))
                .add(String.valueOf(recoveredVariant.prevGCContent()))
                .add(String.valueOf(recoveredVariant.previous()))

                .add(String.valueOf(recoveredVariant.nextGCContent()))
                .add(String.valueOf(recoveredVariant.next()))
                .add(String.valueOf(recoveredVariant.variant()))
                .add(String.valueOf(recoveredVariant.orientation()))
                .add(String.valueOf(recoveredVariant.qual()))
                .add(String.valueOf(recoveredVariant.filter()))
                .add(String.valueOf(recoveredVariant.mate()))
                .add(String.valueOf(recoveredVariant.mateOrientation()))
                .add(String.valueOf(recoveredVariant.mateMinStart()))
                .add(String.valueOf(recoveredVariant.mateMaxStart()))
                .add(String.valueOf(recoveredVariant.mateSupport()))
                .toString();
    }
}

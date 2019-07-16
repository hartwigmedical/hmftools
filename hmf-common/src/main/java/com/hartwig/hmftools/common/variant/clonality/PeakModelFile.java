package com.hartwig.hmftools.common.variant.clonality;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class PeakModelFile {

    ;
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private static final String DELIMITER = "\t";
    private static final String EXTENSION = ".purple.somatic.clonality.tsv";

    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull final String filename, @NotNull final List<PeakModel> bafs) throws IOException {
        Files.write(new File(filename).toPath(), toLines(bafs));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<PeakModel> purity) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(PeakModelFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("peak")
                .add("bucket")
                .add("bucketWeight")
                .add("peakAvgWeight")
                .add("isValid")
                .add("isSubclonal")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final PeakModel ratio) {
        return new StringJoiner(DELIMITER).add(FORMAT.format(ratio.peak()))
                .add(FORMAT.format(ratio.bucket()))
                .add(FORMAT.format(ratio.bucketWeight()))
                .add(FORMAT.format(ratio.peakAvgWeight()))
                .add(String.valueOf(ratio.isValid()))
                .add(String.valueOf(ratio.isSubclonal()))
                .toString();
    }

}

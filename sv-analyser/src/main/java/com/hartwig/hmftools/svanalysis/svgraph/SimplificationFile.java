package com.hartwig.hmftools.svanalysis.svgraph;

import com.google.common.collect.Lists;
import org.jetbrains.annotations.NotNull;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.StringJoiner;

public enum SimplificationFile {
    ;
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00");
    private static final String DELIMITER = "\t";
    static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".simplification.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + EXTENSION;
    }

    public static void write(@NotNull final String filename, @NotNull Collection<Simplification> simplifications) throws IOException {
        Files.write(new File(filename).toPath(), toLines(simplifications));
    }

    @NotNull
    static List<String> toLines(@NotNull final Collection<Simplification> variants) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        int id = 0;
        for (Simplification s : variants) {
            lines.addAll(toFlatStrings(s, id++));
        }
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "")
                .add("uid")
                .add("type")
                .add("ploidy")
                .add("svid")
                .add("anchorcn")
                .add("refcn")
                .add("svploidy")
                .add("deltacn")
                .add("deltasv")
                .add("othersvploidy")
                .toString();
    }

    @NotNull
    private static List<String> toFlatStrings(@NotNull final Simplification simplification, int id) {
        List<String> list = new ArrayList<>();
        for (BreakendConsistency bc : simplification.consistency()) {
            list.add(new StringJoiner(DELIMITER)
                    .add(Integer.toString(id))
                    .add(String.valueOf(simplification.type()))
                    .add(FORMAT.format(simplification.ploidy()))
                    .add(Integer.toString(bc.sv().primaryKey()))
                    .add(FORMAT.format(bc.anchorPloidy()))
                    .add(FORMAT.format(bc.referencePathPloidy()))
                    .add(FORMAT.format(bc.sv().ploidy()))
                    .add(FORMAT.format(bc.copyNumberDelta()))
                    .add(FORMAT.format(bc.eventDelta()))
                    .add(FORMAT.format(bc.otherSvPloidy()))
                    .toString());
        }
        return list;
    }
}

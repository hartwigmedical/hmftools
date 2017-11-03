package com.hartwig.hmftools.common.amber;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public enum AmberBAFFile {
    ;
    private static final DecimalFormat FORMAT = new DecimalFormat("0.0000");

    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "Chr";
    private static final String PURPLE_EXTENSION = ".purple.baf";
    private static final String AMBER_EXTENSION = ".amber.baf";



    @Deprecated
    public static String generatePurpleFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + PURPLE_EXTENSION;
    }

    public static String generateAmberFilename(@NotNull final String basePath, @NotNull final String sample) {
        return basePath + File.separator + sample + AMBER_EXTENSION;
    }

    @NotNull
    public static Multimap<String, AmberBAF> read(@NotNull final String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull Multimap<String, AmberBAF> bafs) throws IOException {
        List<AmberBAF> sortedBafs = Lists.newArrayList(bafs.values());
        Collections.sort(sortedBafs);
        write(filename, sortedBafs);
    }

    public static void write(@NotNull final String filename, @NotNull List<AmberBAF> bafs) throws IOException {
        Files.write(new File(filename).toPath(), toLines(bafs));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<AmberBAF> purity) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(AmberBAFFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, "", "").add("Chromosome")
                .add("Position")
                .add("TumorBAF")
                .add("TumorModifiedBAF")
                .add("TumorDepth")
                .add("NormalBAF")
                .add("NormalModifiedBAF")
                .add("NormalDepth")
                .toString();
    }

    @NotNull
    private static String toString(@NotNull final AmberBAF ratio) {
        return new StringJoiner(DELIMITER).add(String.valueOf(ratio.chromosome()))
                .add(String.valueOf(ratio.position()))
                .add(FORMAT.format(ratio.tumorBAF()))
                .add(FORMAT.format(ratio.tumorModifiedBAF()))
                .add(String.valueOf(ratio.tumorDepth()))
                .add(FORMAT.format(ratio.normalBAF()))
                .add(FORMAT.format(ratio.normalModifiedBAF()))
                .add(String.valueOf(ratio.normalDepth()))
                .toString();
    }

    @NotNull
    private static Multimap<String, AmberBAF> fromLines(@NotNull List<String> lines) {
        Multimap<String, AmberBAF> result = ArrayListMultimap.create();
        for (String line : lines) {
            if (!line.startsWith(HEADER_PREFIX)) {
                final AmberBAF region = fromString(line);
                result.put(region.chromosome(), region);
            }
        }
        return result;
    }

    @NotNull
    private static AmberBAF fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        ImmutableAmberBAF.Builder builder = ImmutableAmberBAF.builder()
                .chromosome(values[0])
                .position(Long.valueOf(values[1]))
                .tumorBAF(Double.valueOf(values[2]))
                .tumorDepth(0)
                .normalBAF(0.5)
                .normalDepth(0);

        if (values.length == 8) {
            builder.tumorDepth(Integer.valueOf(values[4])).normalBAF(Double.valueOf(values[5])).normalDepth(Integer.valueOf(values[7]));
        }

        return builder.build();
    }
}

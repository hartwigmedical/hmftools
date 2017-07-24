package com.hartwig.hmftools.common.purple.ratio;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum ReadRatioFile {
    ;

    private static final String DELIMITER = "\t";
    private static final String HEADER_PREFIX = "#";
    private static final String EXTENSION = ".purple.ratio";

    @NotNull
    public static List<ReadRatio> read(@NotNull final String basePath, @NotNull final String sample) throws IOException {
        final String filePath = basePath + File.separator + sample + EXTENSION;
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String basePath, @NotNull final String sample, @NotNull List<ReadRatio> copyNumbers)
            throws IOException {
        final String filePath = basePath + File.separator + sample + EXTENSION;
        Files.write(new File(filePath).toPath(), toLines(copyNumbers));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<ReadRatio> purity) {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        purity.stream().map(ReadRatioFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header() {
        return new StringJoiner(DELIMITER, HEADER_PREFIX, "").add("chromosome").add("position").add("ratio").toString();
    }

    @NotNull
    private static String toString(@NotNull final ReadRatio ratio) {
        return new StringJoiner(DELIMITER).add(String.valueOf(ratio.chromosome()))
                .add(String.valueOf(ratio.position()))
                .add(String.valueOf(ratio.ratio()))
                .toString();
    }

    @NotNull
    private static List<ReadRatio> fromLines(@NotNull List<String> lines) {
        return lines.stream().filter(x -> !x.startsWith(HEADER_PREFIX)).map(ReadRatioFile::fromString).collect(toList());
    }

    @NotNull
    private static ReadRatio fromString(@NotNull final String line) {
        String[] values = line.split(DELIMITER);
        return ImmutableReadRatio.builder()
                .chromosome(values[0])
                .position(Long.valueOf(values[1]))
                .ratio(Double.valueOf(values[2]))
                .build();
    }
}

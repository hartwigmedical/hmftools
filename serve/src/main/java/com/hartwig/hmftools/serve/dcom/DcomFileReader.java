package com.hartwig.hmftools.serve.dcom;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.jetbrains.annotations.NotNull;

public final class DcomFileReader {

    private static final String DELIMITER = "\t";

    private DcomFileReader() {
    }

    @NotNull
    public static List<DcomEntry> readDcomFile(@NotNull String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static List<DcomEntry> fromLines(@NotNull List<String> lines) {
        return lines.stream().skip(1).map(DcomFileReader::fromLine).collect(toList());
    }

    @NotNull
    private static DcomEntry fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        String hgvs = values[0];
        String transcript = hgvs.substring(0, hgvs.indexOf(":"));
        return ImmutableDcomEntry.builder().transcript(transcript).gene(values[7]).proteinImpact(values[9]).build();
    }
}

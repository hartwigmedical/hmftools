package com.hartwig.hmftools.serve.sources.hartwig.cohort;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.serve.sources.hartwig.HartwigProteinInterpreter;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class HartwigCohortFileReader {

    private static final String DELIMITER = "\t";

    private HartwigCohortFileReader() {
    }

    @NotNull
    public static List<HartwigCohortEntry> readCohortFile(@NotNull String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static List<HartwigCohortEntry> fromLines(@NotNull List<String> lines) {
        return lines.stream().skip(1).map(HartwigCohortFileReader::fromLine).collect(toList());
    }

    @NotNull
    private static HartwigCohortEntry fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        String proteinAnnotation = values.length > 6 ? values[6] : Strings.EMPTY;

        return ImmutableHartwigCohortEntry.builder()
                .chromosome(values[0])
                .position(Long.parseLong(values[1]))
                .ref(values[2])
                .alt(values[3])
                .gene(values[4])
                .transcript(values[5])
                .proteinAnnotation(HartwigProteinInterpreter.interpretProtein(proteinAnnotation))
                .build();
    }
}

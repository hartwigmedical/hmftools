package com.hartwig.hmftools.serve.sources.hartwig.curated;

import static java.util.stream.Collectors.toList;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.serve.sources.hartwig.HartwigProteinInterpreter;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class HartwigCuratedFileReader {

    private static final String DELIMITER = "\t";

    private HartwigCuratedFileReader() {
    }

    @NotNull
    public static List<HartwigCuratedEntry> readCuratedFile(@NotNull String fileName) throws IOException {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static List<HartwigCuratedEntry> fromLines(@NotNull List<String> lines) {
        return lines.stream().skip(1).map(HartwigCuratedFileReader::fromLine).collect(toList());
    }

    @NotNull
    private static HartwigCuratedEntry fromLine(@NotNull String line) {
        String[] values = line.split(DELIMITER);

        String proteinAnnotation = values.length > 6 ? values[6] : Strings.EMPTY;

        return ImmutableHartwigCuratedEntry.builder()
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

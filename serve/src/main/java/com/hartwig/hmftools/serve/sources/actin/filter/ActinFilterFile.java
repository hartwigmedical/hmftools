package com.hartwig.hmftools.serve.sources.actin.filter;

import static com.hartwig.hmftools.serve.actionability.util.ActionableFileFunctions.FIELD_DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class ActinFilterFile {

    private ActinFilterFile() {
    }

    @NotNull
    public static List<ActinFilterEntry> read(@NotNull String actinFilterTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actinFilterTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static List<ActinFilterEntry> fromLines(@NotNull List<String> lines) {
        List<ActinFilterEntry> filterEntries = Lists.newArrayList();
        for (String line : lines) {
            filterEntries.add(fromLine(line));
        }
        return filterEntries;
    }

    @NotNull
    private static ActinFilterEntry fromLine(@NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableActinFilterEntry.builder()
                .type(ActinFilterType.valueOf(values[0]))
                .value(values[1])
                .build();
    }
}

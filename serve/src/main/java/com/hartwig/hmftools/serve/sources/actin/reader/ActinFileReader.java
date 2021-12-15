package com.hartwig.hmftools.serve.sources.actin.reader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class ActinFileReader {

    private static final String MAIN_FIELD_DELIMITER = "\t";
    private static final String PARAMETER_DELIMITER = ",";

    private ActinFileReader() {
    }

    public static List<ActinEntry> read(@NotNull String actinTrialTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(actinTrialTsv).toPath()));
    }

    @NotNull
    @VisibleForTesting
    static List<ActinEntry> fromLines(@NotNull List<String> lines) {
        List<ActinEntry> trials = Lists.newArrayList();

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            trials.add(fromString(line));
        }

        return trials;
    }

    @NotNull
    private static ActinEntry fromString(@NotNull String line) {
        String[] values = line.split(MAIN_FIELD_DELIMITER);

        return ImmutableActinEntry.builder()
                .trial(values[0])
                .rule(ActinRule.valueOf(values[1]))
                .parameters(toParameters(values[2]))
                .build();
    }

    @NotNull
    private static List<String> toParameters(@NotNull String parameterString) {
        if (parameterString.isEmpty()) {
            return Lists.newArrayList();
        }

        List<String> parameters = Lists.newArrayList();
        for (String parameter : parameterString.split(PARAMETER_DELIMITER)) {
            parameters.add(parameter.trim());
        }
        return parameters;
    }
}

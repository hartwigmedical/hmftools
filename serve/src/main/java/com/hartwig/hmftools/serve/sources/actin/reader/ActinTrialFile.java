package com.hartwig.hmftools.serve.sources.actin.reader;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.serve.sources.actin.ActinTrial;
import com.hartwig.hmftools.serve.sources.actin.ImmutableActinTrial;

import org.jetbrains.annotations.NotNull;

public class ActinTrialFile {

    public static final String MAIN_FIELD_DELIMITER = "\t";

    private ActinTrialFile() {
    }

    public static List<ActinTrial> read(@NotNull String actinTrialsTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(actinTrialsTsv).toPath()));
    }

    @NotNull
    @VisibleForTesting
    static List<ActinTrial> fromLines(@NotNull List<String> lines) {
        List<ActinTrial> trials = Lists.newArrayList();

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            trials.add(fromString(line));
        }

        return trials;
    }

    @NotNull
    private static ActinTrial fromString(@NotNull String line) {
        String[] values = line.split(MAIN_FIELD_DELIMITER);

        return ImmutableActinTrial.builder()
                .trialId(values[0])
                .cohortId(values[1])
                .rule(values[2])
                .parameters(values[3])
                .build();
    }

}

package com.hartwig.hmftools.serve.sources.actin.reader;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ActinFileReader {

    private static final String MAIN_FIELD_DELIMITER = "\t";

    private ActinFileReader() {
    }

    @NotNull
    public static List<ActinEntry> read(@NotNull String actinTrialTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(actinTrialTsv).toPath()));
    }

    @NotNull
    private static List<ActinEntry> fromLines(@NotNull List<String> lines) {
        List<ActinEntry> trials = Lists.newArrayList();

        Map<String, Integer> fields = createFieldsIndexMap(lines.get(0), MAIN_FIELD_DELIMITER);
        for (String line : lines.subList(1, lines.size())) {
            trials.add(fromString(fields, line));
        }

        return trials;
    }

    @NotNull
    private static ActinEntry fromString(@NotNull Map<String, Integer> fields, @NotNull String line) {
        String[] values = line.split(MAIN_FIELD_DELIMITER, -1);

        return ImmutableActinEntry.builder()
                .trial(values[fields.get("trial")])
                .rule(ActinRule.valueOf(values[fields.get("rule")]))
                .gene(emptyToNull(values[fields.get("gene")]))
                .mutation(emptyToNull(values[fields.get("mutation")]))
                .build();
    }

    @Nullable
    private static String emptyToNull(@NotNull String value) {
        return !value.isEmpty() ? value : null;
    }
}

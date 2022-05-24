package com.hartwig.hmftools.rose.actionability;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.FileReaderUtils;

import org.jetbrains.annotations.NotNull;

public final class ActionabilityFileReader {

    private static final String MAIN_FIELD_DELIMITER = "\t";

    private ActionabilityFileReader() {
    }

    @NotNull
    public static List<ActionabilityEntry> read(@NotNull String actionabilityTsv) throws IOException {
        return fromLines(Files.readAllLines(new File(actionabilityTsv).toPath()));
    }

    @NotNull
    private static List<ActionabilityEntry> fromLines(@NotNull List<String> lines) {
        List<ActionabilityEntry> trials = Lists.newArrayList();

        Map<String, Integer> fields = FileReaderUtils.createFieldsIndexMap(lines.get(0), MAIN_FIELD_DELIMITER);
        for (String line : lines.subList(1, lines.size())) {
            trials.add(fromString(fields, line));
        }
        return trials;
    }

    @NotNull
    private static ActionabilityEntry fromString(@NotNull Map<String, Integer> fields, @NotNull String line) {
        String[] values = line.split(MAIN_FIELD_DELIMITER, -1);

        return ImmutableActionabilityEntry.builder()
                .match(values[fields.get("match")])
                .type(TypeAlteration.toType(values[fields.get("type_alteration")]))
                .condition(Condition.toCondition(values[fields.get("condition")]))
                .conclusion(values[fields.get("conclusion")])
                .build();
    }
}
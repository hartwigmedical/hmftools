package com.hartwig.hmftools.summon.actionability;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.FileReaderUtils;
import com.hartwig.hmftools.summon.SummonApplication;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class ActionabilityFileReader {
    private static final Logger LOGGER = LogManager.getLogger(ActionabilityFileReader.class);

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
                .gene(values[fields.get("Gene")])
                .type(Type.valueOf(values[fields.get("aberration type")].toUpperCase().replace(" ", "_")))
                .conclusion(values[fields.get("conclusion")])
                .build();
    }
}
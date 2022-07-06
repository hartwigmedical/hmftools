package com.hartwig.hmftools.serve.curation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class DrugClassFactory {

    private static final String FIELD_DELIMITER = "\t";

    @NotNull
    public static Map<DrugClassKey, DrugClasses> read(@NotNull String drugClassCurationTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(drugClassCurationTsv).toPath());
        // Skip header
        return fromLines(lines.subList(1, lines.size()));
    }

    @NotNull
    private static Map<DrugClassKey, DrugClasses> fromLines(@NotNull List<String> lines) {
        Map<DrugClassKey, DrugClasses> curatedEntries = Maps.newHashMap();

        for (String line : lines) {
            String[] values = line.split(FIELD_DELIMITER);
            DrugClassKey drugClassKey =
                    ImmutableDrugClassKey.builder().treatment(values[0]).drugClass(values[1]).matchEvent(values[2]).build();
            curatedEntries.put(drugClassKey, ImmutableDrugClasses.builder().drugClassKey(drugClassKey).curatedDrugClass(values[3]).build());
        }
        return curatedEntries;
    }
}
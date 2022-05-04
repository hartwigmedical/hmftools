package com.hartwig.hmftools.serve.sources.actin.filter;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class ActinFilterFile {

    private static final String FIELD_DELIMITER = "\t";

    private ActinFilterFile() {
    }

    @NotNull
    public static List<ActinFilterEntry> read(@NotNull String actinFilterTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(actinFilterTsv).toPath());

        Map<String, Integer> fields = createFieldsIndexMap(lines.get(0), FIELD_DELIMITER);
        List<ActinFilterEntry> filterEntries = Lists.newArrayList();
        for (String line : lines.subList(1, lines.size())) {
            filterEntries.add(fromLine(fields, line));
        }
        return filterEntries;
    }

    @NotNull
    private static ActinFilterEntry fromLine(@NotNull Map<String, Integer> fields, @NotNull String line) {
        String[] values = line.split(FIELD_DELIMITER);

        return ImmutableActinFilterEntry.builder()
                .type(ActinFilterType.valueOf(values[fields.get("filterType")]))
                .value(values[fields.get("value")])
                .build();
    }
}

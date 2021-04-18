package com.hartwig.hmftools.serve.curation;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public final class DoidLookupFactory {

    private static final String FIELD_DELIMITER = "\t";
    private static final String DOID_DELIMITER = ";";

    private DoidLookupFactory() {
    }

    @NotNull
    public static DoidLookup buildFromMappingTsv(@NotNull String configTsvPath) throws IOException {
        Map<String, Set<String>> cancerTypeToDoidsMapping = Maps.newHashMap();

        List<String> lines = Files.readAllLines(new File(configTsvPath).toPath());
        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_DELIMITER);
            cancerTypeToDoidsMapping.put(parts[0], Sets.newHashSet(parts[1].split(DOID_DELIMITER)));
        }

        return new DoidLookup(cancerTypeToDoidsMapping);
    }
}

package com.hartwig.hmftools.common.actionability.cancertype;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class CancerTypeToDOIDMapper {

    private static final Logger LOGGER = LogManager.getLogger(CancerTypeToDOIDMapper.class);

    private static final String LINE_DELIMITER = "\t";
    private static final String DOID_DELIMITER = ";";

    @NotNull
    private final Map<String, Set<String>> doidsPerCancerType;

    @NotNull
    public static CancerTypeToDOIDMapper createFromFile(@NotNull String knowledgebaseCancerTypesPath) throws IOException {
        final Map<String, Set<String>> cancerTypeMappings = Maps.newHashMap();
        final List<String> cancerTypeMappingLines = Files.readAllLines(new File(knowledgebaseCancerTypesPath).toPath());

        // KODU: Skip header line
        for (String cancerTypeMappingLine : cancerTypeMappingLines.subList(1, cancerTypeMappingLines.size())) {
            CancerTypeToDOIDMappingEntry entry = fromLine(cancerTypeMappingLine);
            cancerTypeMappings.put(entry.cancerType(), entry.doids());
        }

        return new CancerTypeToDOIDMapper(cancerTypeMappings);
    }

    @VisibleForTesting
    CancerTypeToDOIDMapper(@NotNull final Map<String, Set<String>> doidsPerCancerType) {
        this.doidsPerCancerType = doidsPerCancerType;
    }

    @NotNull
    private static CancerTypeToDOIDMappingEntry fromLine(@NotNull String line) {
        String[] values = line.split(LINE_DELIMITER);
        Set<String> doids = values.length > 1 ? toDOIDSet(values[1].trim()) : Sets.newHashSet();
        return ImmutableCancerTypeToDOIDMappingEntry.builder().cancerType(values[0].trim()).doids(doids).build();
    }

    @NotNull
    private static Set<String> toDOIDSet(@NotNull String doidValue) {
        String[] values = doidValue.split(DOID_DELIMITER);
        Set<String> doids = Sets.newHashSet();
        for (String value : values) {
            // KODU: Expected format is "Doid(value=NNN)"
            String doid = value.substring(value.indexOf("=") + 1, value.indexOf(")"));
            doids.add(doid.trim());
        }

        return doids;
    }

    @Nullable
    public Set<String> findDoids(@NotNull String cancerType) {
        Set<String> doids = doidsPerCancerType.get(cancerType);
        if (doids == null) {
            LOGGER.warn("Could not resolve cancer type in DOID mapping: " + cancerType);
        }
        return doids;
    }
}

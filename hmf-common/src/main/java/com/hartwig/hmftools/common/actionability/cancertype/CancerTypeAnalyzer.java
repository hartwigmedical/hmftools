package com.hartwig.hmftools.common.actionability.cancertype;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CancerTypeAnalyzer {
    private static final Logger LOGGER = LogManager.getLogger(CancerTypeAnalyzer.class);

    private static final String LINE_DELIMITER = "\t";
    private static final String DOID_DELIMITER = ";";

    @NotNull
    private final List<CancerTypeToDOIDMappingEntry> cancerTypeMappings;
    @NotNull
    private final PrimaryTumorToDOIDMapping primaryTumorToDOIDMapping;

    public CancerTypeAnalyzer(@NotNull final List<CancerTypeToDOIDMappingEntry> cancerTypeMappings,
            @NotNull final PrimaryTumorToDOIDMapping primaryTumorToDOIDMapping) {
        this.cancerTypeMappings = cancerTypeMappings;
        this.primaryTumorToDOIDMapping = primaryTumorToDOIDMapping;
    }

    @NotNull
    public static CancerTypeAnalyzer loadFromFile(@NotNull String knowledgebaseCancerTypesPath) throws IOException {
        final List<CancerTypeToDOIDMappingEntry> cancerTypeMappings = Lists.newArrayList();
        final List<String> cancerTypeMappingLines = Files.readAllLines(new File(knowledgebaseCancerTypesPath).toPath());

        // KODU: Skip header line
        for (String cancerTypeMappingLine : cancerTypeMappingLines.subList(1, cancerTypeMappingLines.size())) {
            cancerTypeMappings.add(fromLine(cancerTypeMappingLine));
        }

        return new CancerTypeAnalyzer(cancerTypeMappings, PrimaryTumorToDOIDMapping.createFromResource());
    }

    @NotNull
    private static CancerTypeToDOIDMappingEntry fromLine(@NotNull String line) {
        String[] values = line.split(LINE_DELIMITER);
        return ImmutableCancerTypeToDOIDMappingEntry.builder().cancerType(values[0].trim()).doids(toDOIDSet(values[1].trim())).build();
    }

    @NotNull
    private static Set<Integer> toDOIDSet(@NotNull String doidValue) {
        String[] values = doidValue.split(DOID_DELIMITER);
        Set<Integer> doids = Sets.newHashSet();
        for (String value : values) {
            // KODU: Expected format is "Doid(value=NNN)"
            String doid = value.substring(value.indexOf("=") + 1, value.indexOf(")"));
            doids.add(Integer.valueOf(doid));
        }

        return doids;
    }

    public boolean isCancerTypeMatch(@NotNull String knowledgebaseCancerType, @Nullable String primaryTumorLocation) {
        if (primaryTumorLocation == null) {
            return false;
        }

        Set<Integer> doidsForPrimaryTumorLocation = primaryTumorToDOIDMapping.findDoids(primaryTumorLocation);
        if (doidsForPrimaryTumorLocation == null) {
            return false;
        }

        Set<Integer> doidsForCancerType = findDoidsForCancerType(knowledgebaseCancerType, cancerTypeMappings);
        return doidsForCancerType != null && !Collections.disjoint(doidsForPrimaryTumorLocation, doidsForCancerType);
    }

    @Nullable
    private static Set<Integer> findDoidsForCancerType(@NotNull String cancerType,
            @NotNull List<CancerTypeToDOIDMappingEntry> cancerTypeMappings) {
        for (CancerTypeToDOIDMappingEntry mapping : cancerTypeMappings) {
            if (mapping.cancerType().equalsIgnoreCase(cancerType)) {
                return mapping.doids();
            }
        }

        LOGGER.warn("Could not resolve cancer type in DOID mapping: " + cancerType);
        return null;
    }
}

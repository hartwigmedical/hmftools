package com.hartwig.hmftools.patientdb.clinical.curators;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableCuratedPrimaryTumor;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PrimaryTumorCurator implements CleanableCurator {

    private static final Logger LOGGER = LogManager.getLogger(PrimaryTumorCurator.class);

    private static final String FIELD_DELIMITER = "\t";
    private static final String STRING_DELIMITER = ";";

    @NotNull
    private final Map<String, CuratedPrimaryTumor> primaryTumorMap;
    @NotNull
    private final Map<String, String> tumorLocationOverridesMap;
    @NotNull
    private final Set<String> unusedSearchTerms;
    @NotNull
    private final DoidNodesResolver doidNodesResolver;

    public PrimaryTumorCurator(@NotNull String tumorLocationMappingTsv, @NotNull String tumorLocationOverridesTsv,
            @NotNull DoidNodesResolver doidNodesResolver) throws IOException {
        tumorLocationOverridesMap = loadFromOverridesTsv(tumorLocationOverridesTsv);
        this.doidNodesResolver = doidNodesResolver;
        primaryTumorMap = loadFromMappingTsv(tumorLocationMappingTsv);


        // Need to create a copy of the key set so that we can remove elements from it without affecting the curation.
        unusedSearchTerms = Sets.newHashSet(primaryTumorMap.keySet());
    }

    @NotNull
    private Map<String, CuratedPrimaryTumor> loadFromMappingTsv(@NotNull String tumorLocationMappingTsv) throws IOException {
        Map<String, CuratedPrimaryTumor> primaryTumorMap = Maps.newHashMap();

        List<String> lines = Files.readAllLines(new File(tumorLocationMappingTsv).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_DELIMITER);
            String searchTerm = parts[0];
            String location = parts[1];
            String subLocation = parts[2];
            String type = parts[3];
            String subType = parts[4];
            String extraDetails = parts[5];
            List<String> doids = Lists.newArrayList(parts[6].split(STRING_DELIMITER));
            List<String> snomedConceptIds = Lists.newArrayList(parts[7].split(STRING_DELIMITER));
            List<DoidNode> resolvedDoidNodes = doidNodesResolver.resolveDoidNodes(doids);

            primaryTumorMap.put(searchTerm,
                    ImmutableCuratedPrimaryTumor.builder()
                            .searchTerm(searchTerm)
                            .isOverridden(false)
                            .location(location)
                            .subLocation(subLocation)
                            .type(type)
                            .subType(subType)
                            .extraDetails(extraDetails)
                            .doidNodes(resolvedDoidNodes)
                            .snomedConceptIds(snomedConceptIds)
                            .build());
        }
        return primaryTumorMap;
    }

    @NotNull
    private static Map<String, String> loadFromOverridesTsv(@NotNull String tumorLocationOverridesTsv) throws IOException {
        Map<String, String> tumorLocationOverridesMap = Maps.newHashMap();
        List<String> lines = Files.readAllLines(new File(tumorLocationOverridesTsv).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_DELIMITER);
            tumorLocationOverridesMap.put(parts[0], parts[1]);
        }
        return tumorLocationOverridesMap;
    }

    @NotNull
    public CuratedPrimaryTumor search(@NotNull String patientIdentifier, @Nullable String searchTerm) {
        String effectiveSearchTerm = searchTerm;
        String override = tumorLocationOverridesMap.get(patientIdentifier);
        if (override != null) {
            LOGGER.info("  Overriding tumor location for {} from '{}' to '{}'", patientIdentifier, searchTerm, override);
            effectiveSearchTerm = override;
        }

        if (effectiveSearchTerm != null && !effectiveSearchTerm.isEmpty()) {
            unusedSearchTerms.remove(effectiveSearchTerm);
            CuratedPrimaryTumor result = primaryTumorMap.get(effectiveSearchTerm);

            if (result != null) {
                return ImmutableCuratedPrimaryTumor.builder().from(result).isOverridden(override != null).build();
            }
        }

        return ImmutableCuratedPrimaryTumor.builder().isOverridden(false).searchTerm(effectiveSearchTerm).build();
    }

    @NotNull
    @Override
    public Set<String> unusedSearchTerms() {
        return unusedSearchTerms;
    }
}

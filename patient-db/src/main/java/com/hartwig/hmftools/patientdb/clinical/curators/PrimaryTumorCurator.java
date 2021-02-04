package com.hartwig.hmftools.patientdb.clinical.curators;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.clinical.data.CuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.clinical.data.ImmutableCuratedPrimaryTumor;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PrimaryTumorCurator implements CleanableCurator {

    private static final String FIELD_DELIMITER = "\t";
    private static final String STRING_DELIMITER = ";";

    @NotNull
    private final Map<String, CuratedPrimaryTumor> primaryTumorMap = Maps.newHashMap();
    @NotNull
    private final Set<String> unusedSearchTerms;

    public PrimaryTumorCurator(@NotNull String tumorLocationMappingTsv, @NotNull List<DoidNode> doidNodes) throws IOException {
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

            primaryTumorMap.put(searchTerm,
                    ImmutableCuratedPrimaryTumor.builder()
                            .searchTerm(searchTerm)
                            .location(location)
                            .subLocation(subLocation)
                            .type(type)
                            .subType(subType)
                            .extraDetails(extraDetails)
                            .doidNodes(resolveDoidNodes(doidNodes, doids))
                            .snomedConceptIds(snomedConceptIds)
                            .build());
        }

        // Need to create a copy of the key set so that we can remove elements from it without affecting the curation.
        unusedSearchTerms = Sets.newHashSet(primaryTumorMap.keySet());
    }

    @NotNull
    public CuratedPrimaryTumor search(@Nullable String searchTerm) {
        if (searchTerm != null && !searchTerm.isEmpty()) {
            unusedSearchTerms.remove(searchTerm);
            CuratedPrimaryTumor result = primaryTumorMap.get(searchTerm);

            if (result != null) {
                return result;
            }
        }

        return ImmutableCuratedPrimaryTumor.builder().searchTerm(searchTerm).build();
    }

    @NotNull
    @Override
    public Set<String> unusedSearchTerms() {
        return unusedSearchTerms;
    }

    @NotNull
    @VisibleForTesting
    static List<DoidNode> resolveDoidNodes(@NotNull List<DoidNode> doidNodes, @NotNull List<String> doidsToResolve) {
        List<DoidNode> resolvedDoidNodes = Lists.newArrayList();
        for (String doid : doidsToResolve) {
            for (DoidNode doidNode : doidNodes) {
                if (doidNode.doid().equals(doid)) {
                    resolvedDoidNodes.add(doidNode);
                }
            }
        }
        return resolvedDoidNodes;
    }
}

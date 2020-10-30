package com.hartwig.hmftools.patientdb.curators;

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
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorLocationCurator implements CleanableCurator {
    private static final Logger LOGGER = LogManager.getLogger(TumorLocationCurator.class);

    private static final String FIELD_DELIMITER = "\t";
    private static final String DOID_DELIMITER = ";";

    @NotNull
    private final Map<String, CuratedTumorLocation> tumorLocationMap = Maps.newHashMap();
    @NotNull
    private final Set<String> unusedSearchTerms;

    public TumorLocationCurator(@NotNull String tumorLocationMappingTsv, @NotNull List<DoidNode> doidNodes) throws IOException {
        List<String> lines = Files.readAllLines(new File(tumorLocationMappingTsv).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_DELIMITER);
            String searchTerm = parts[0];
            String primaryTumorLocation = parts[1];
            String primaryTumorSubLocation = parts.length > 2 ? parts[2] : Strings.EMPTY;
            String primaryTumorType = parts.length > 3 ? parts[3] : Strings.EMPTY;
            String primaryTumorSubType = parts.length > 4 ? parts[4] : Strings.EMPTY;
            String primaryTumorExtraDetails = parts.length > 5 ? parts[5] : Strings.EMPTY;
            List<String> doids = parts.length > 6 ? Lists.newArrayList(parts[6].split(DOID_DELIMITER)) : Lists.newArrayList();

            tumorLocationMap.put(searchTerm,
                    ImmutableCuratedTumorLocation.builder()
                            .primaryTumorLocation(primaryTumorLocation)
                            .primaryTumorSubLocation(primaryTumorSubLocation)
                            .primaryTumorType(primaryTumorType)
                            .primaryTumorSubType(primaryTumorSubType)
                            .primaryTumorExtraDetails(primaryTumorExtraDetails)
                            .doidNodes(resolveDoidNodes(doidNodes, doids))
                            .searchTerm(searchTerm)
                            .build());
        }

        // Need to create a copy of the key set so that we can remove elements from it without affecting the curation.
        unusedSearchTerms = Sets.newHashSet(tumorLocationMap.keySet());
    }

    @NotNull
    @VisibleForTesting
    public static List<DoidNode> resolveDoidNodes(@NotNull List<DoidNode> doidNodes, @NotNull List<String> doidNodesToResolve) {
        LOGGER.info(doidNodes);
        List<DoidNode> resolvedDoidNodes = Lists.newArrayList();
        for (String doid : doidNodesToResolve) {
            for (DoidNode doidNode : doidNodes) {
                if (doidNode.doid().equals(doid)) {
                    resolvedDoidNodes.add(doidNode);
                }
            }
        }
        return resolvedDoidNodes;

    }

    @NotNull
    public CuratedTumorLocation search(@Nullable String searchTerm) {
        if (searchTerm != null && !searchTerm.isEmpty()) {
            unusedSearchTerms.remove(searchTerm);
            CuratedTumorLocation result = tumorLocationMap.get(searchTerm);

            if (result != null) {
                return result;
            }
        }

        return ImmutableCuratedTumorLocation.builder().searchTerm(searchTerm).build();
    }

    @NotNull
    @Override
    public Set<String> unusedSearchTerms() {
        return unusedSearchTerms;
    }
}

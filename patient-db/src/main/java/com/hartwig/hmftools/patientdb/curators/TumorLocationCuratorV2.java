package com.hartwig.hmftools.patientdb.curators;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.doid.DoidEntry;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocationV2;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocationV2;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorLocationCuratorV2 implements CleanableCurator {

    private static final String FIELD_DELIMITER = "\t";
    private static final String DOID_DELIMITER = ";";

    @NotNull
    private final Map<String, CuratedTumorLocationV2> tumorLocationMap = Maps.newHashMap();
    @NotNull
    private final Set<String> unusedSearchTerms;

    public TumorLocationCuratorV2(@NotNull String tumorLocationV2MappingTsv, @NotNull List<DoidEntry> doidEntries) throws IOException {
        List<String> lines = Files.readAllLines(new File(tumorLocationV2MappingTsv).toPath());

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
                    ImmutableCuratedTumorLocationV2.builder()
                            .primaryTumorLocation(primaryTumorLocation)
                            .primaryTumorSubLocation(primaryTumorSubLocation)
                            .primaryTumorType(primaryTumorType)
                            .primaryTumorSubType(primaryTumorSubType)
                            .primaryTumorExtraDetails(primaryTumorExtraDetails)
                            .doidEntries(resolveDoidEntries(doidEntries, doids))
                            .searchTerm(searchTerm)
                            .build());
        }

        // Need to create a copy of the key set so that we can remove elements from it without affecting the curation.
        unusedSearchTerms = Sets.newHashSet(tumorLocationMap.keySet());
    }

    @NotNull
    private static List<DoidEntry> resolveDoidEntries(@NotNull List<DoidEntry> doidEntries, @NotNull List<String> doidsToResolve) {
        List<DoidEntry> resolvedDoidEntries = Lists.newArrayList();
        for (String doid : doidsToResolve) {
            for (DoidEntry doidEntry : doidEntries) {
                if (doidEntry.doid().equals(doid)) {
                    resolvedDoidEntries.add(doidEntry);
                }
            }
        }
        return resolvedDoidEntries;

    }

    @NotNull
    public CuratedTumorLocationV2 search(@Nullable String searchTerm) {
        if (searchTerm != null && !searchTerm.isEmpty()) {
            unusedSearchTerms.remove(searchTerm);
            CuratedTumorLocationV2 result = tumorLocationMap.get(searchTerm);

            if (result != null) {
                return result;
            }
        }

        return ImmutableCuratedTumorLocationV2.builder().searchTerm(searchTerm).build();
    }

    @NotNull
    @Override
    public Set<String> unusedSearchTerms() {
        return unusedSearchTerms;
    }
}

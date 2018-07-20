package com.hartwig.hmftools.patientdb.curators;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocation;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jetbrains.annotations.TestOnly;

public class TumorLocationCurator implements CleanableCurator {

    private static final InputStream TUMOR_LOCATION_MAPPING_RESOURCE =
            TumorLocationCurator.class.getResourceAsStream("/tumor_location_mapping.csv");

    @NotNull
    private final Map<String, CuratedTumorLocation> tumorLocationMap = Maps.newHashMap();
    @NotNull
    private final Set<String> unusedSearchTerms;

    @NotNull
    public static TumorLocationCurator fromProductionResource() throws IOException {
        return new TumorLocationCurator(TUMOR_LOCATION_MAPPING_RESOURCE);
    }

    @VisibleForTesting
    public TumorLocationCurator(@NotNull final InputStream mappingInputStream) throws IOException {
        final CSVParser parser = CSVParser.parse(mappingInputStream, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        for (final CSVRecord record : parser) {
            final String searchTerm = record.get("searchTerm");
            final String primaryTumorLocation = record.get("primaryTumorLocation");
            final String subType = record.get("subType");
            tumorLocationMap.put(searchTerm.toLowerCase(),
                    ImmutableCuratedTumorLocation.of(Utils.capitalize(primaryTumorLocation), Utils.capitalize(subType), searchTerm));
        }
        // KODU: Need to create a copy of the key set so that we can remove elements from it without affecting the curation.
        unusedSearchTerms = Sets.newHashSet(tumorLocationMap.keySet());
    }

    @NotNull
    public CuratedTumorLocation search(@Nullable final String searchTerm) {
        if (searchTerm != null) {
            String effectiveSearchTerm = searchTerm.toLowerCase();
            unusedSearchTerms.remove(effectiveSearchTerm);
            final CuratedTumorLocation result = tumorLocationMap.get(effectiveSearchTerm);

            if (result != null) {
                return result;
            }
        }

        return ImmutableCuratedTumorLocation.of(null, null, searchTerm);
    }

    @NotNull
    @Override
    public Set<String> unusedSearchTerms() {
        return unusedSearchTerms;
    }

    @TestOnly
    @NotNull
    public Set<String> primaryTumorLocations() {
        return tumorLocationMap.entrySet()
                .stream()
                .map(Map.Entry::getValue)
                .map(CuratedTumorLocation::primaryTumorLocation)
                .collect(Collectors.toSet());
    }
}

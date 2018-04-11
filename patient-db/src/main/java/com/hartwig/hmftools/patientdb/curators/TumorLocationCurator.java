package com.hartwig.hmftools.patientdb.curators;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientdb.LoadClinicalData;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.CuratedCancerType;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedCancerType;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorLocationCurator implements CleanableCurator {

    private static final InputStream TUMOR_LOCATION_MAPPING_RESOURCE =
            LoadClinicalData.class.getResourceAsStream("/tumor_location_mapping.csv");

    @NotNull
    private final Map<String, CuratedCancerType> tumorLocationMap = Maps.newHashMap();
    @NotNull
    private final Set<String> unusedSearchTerms;

    @NotNull
    public static TumorLocationCurator fromProductionResource() throws IOException {
        return new TumorLocationCurator(TUMOR_LOCATION_MAPPING_RESOURCE);
    }

    @VisibleForTesting
    TumorLocationCurator(@NotNull final InputStream mappingInputStream) throws IOException {
        final CSVParser parser = CSVParser.parse(mappingInputStream, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        for (final CSVRecord record : parser) {
            final String searchTerm = record.get("searchTerm");
            final String primaryTumorLocation = record.get("primaryTumorLocation");
            final String subType = record.get("subType");
            tumorLocationMap.put(searchTerm.toLowerCase(),
                    ImmutableCuratedCancerType.of(Utils.capitalize(primaryTumorLocation), Utils.capitalize(subType), searchTerm));
        }
        // KODU: Need to create a copy of the key set so that we can remove elements from it without affecting the curation.
        unusedSearchTerms = Sets.newHashSet(tumorLocationMap.keySet());
    }

    @NotNull
    public CuratedCancerType search(@Nullable final String searchTerm) {
        if (searchTerm != null) {
            String effectiveSearchTerm = searchTerm.toLowerCase();
            unusedSearchTerms.remove(effectiveSearchTerm);
            final CuratedCancerType result = tumorLocationMap.get(effectiveSearchTerm);

            if (result != null) {
                return result;
            }
        }

        return ImmutableCuratedCancerType.of(null, null, searchTerm);
    }

    @NotNull
    @Override
    public Set<String> unusedSearchTerms() {
        return unusedSearchTerms;
    }
}

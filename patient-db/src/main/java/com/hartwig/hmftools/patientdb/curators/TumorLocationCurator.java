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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorLocationCurator implements CleanableCurator {

    private static final Logger LOGGER = LogManager.getLogger(TumorLocationCurator.class);
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
            final String location = record.get("primaryTumorLocation");
            final String type = record.get("type");
            final String subType = record.get("subType");
            tumorLocationMap.put(location.toLowerCase(),
                    ImmutableCuratedCancerType.of(Utils.capitalize(type), Utils.capitalize(subType), location));
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

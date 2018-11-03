package com.hartwig.hmftools.common.actionability.cancertype;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Map;
import java.util.Set;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class PrimaryTumorToDOIDMapper {

    private static final Logger LOGGER = LogManager.getLogger(PrimaryTumorToDOIDMapper.class);

    private static final InputStream TUMOR_LOCATION_MAPPING_RESOURCE =
            PrimaryTumorToDOIDMapper.class.getResourceAsStream("/actionability/primary_tumor_locations_mapping.csv");

    private static final String DOID_SEPARATOR = ";";
    private static final String NO_DOIDS_PRESENT = "UNMAPPED";

    @NotNull
    private final Map<String, Set<String>> doidsPerPrimaryTumor;

    @NotNull
    public static PrimaryTumorToDOIDMapper createFromResource() throws IOException {
        final CSVParser parser = CSVParser.parse(TUMOR_LOCATION_MAPPING_RESOURCE, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        Map<String, Set<String>> doidsPerPrimaryTumor = Maps.newHashMap();
        for (final CSVRecord record : parser) {
            final String primaryTumorLocation = record.get("primaryTumorLocation");
            final String doids = record.get("doids");

            doidsPerPrimaryTumor.put(primaryTumorLocation, toSet(doids));
        }

        return new PrimaryTumorToDOIDMapper(doidsPerPrimaryTumor);
    }

    @VisibleForTesting
    PrimaryTumorToDOIDMapper(@NotNull final Map<String, Set<String>> doidsPerPrimaryTumor) {
        this.doidsPerPrimaryTumor = doidsPerPrimaryTumor;
    }

    @NotNull
    private static Set<String> toSet(@NotNull String doidsString) {
        if (doidsString.equals(NO_DOIDS_PRESENT)) {
            return Sets.newHashSet();
        }

        Set<String> doids = Sets.newHashSet();
        String[] values = doidsString.split(DOID_SEPARATOR);

        for (String value : values) {
            doids.add(value.trim());
        }
        return doids;
    }

    @Nullable
    public Set<String> findDoids(@NotNull String primaryTumorLocation) {
        if (primaryTumorLocation.trim().isEmpty()) {
            return null;
        }

        Set<String> doids = doidsPerPrimaryTumor.get(primaryTumorLocation);
        if (doids == null) {
            LOGGER.warn("Could not resolve primary tumor location in DOID mapping: " + primaryTumorLocation);
        }
        return doids;
    }
}

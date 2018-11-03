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
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class PrimaryTumorToDOIDMapping {

    private static final InputStream TUMOR_LOCATION_MAPPING_RESOURCE =
            PrimaryTumorToDOIDMapping.class.getResourceAsStream("/actionability/primary_tumor_locations_mapping.csv");

    private static final String DOID_SEPARATOR = ";";
    private static final String NO_DOIDS_PRESENT = "UNMAPPED";

    @NotNull
    private final Map<String, Set<Integer>> doidsPerPrimaryTumor;

    @NotNull
    public static PrimaryTumorToDOIDMapping createFromResource() throws IOException {
        final CSVParser parser = CSVParser.parse(TUMOR_LOCATION_MAPPING_RESOURCE, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        Map<String, Set<Integer>> doidsPerPrimaryTumor = Maps.newHashMap();
        for (final CSVRecord record : parser) {
            final String primaryTumorLocation = record.get("primaryTumorLocation");
            final String doids = record.get("doids");

            doidsPerPrimaryTumor.put(primaryTumorLocation, toIntegerSet(doids));
        }

        return new PrimaryTumorToDOIDMapping(doidsPerPrimaryTumor);
    }

    @VisibleForTesting
    PrimaryTumorToDOIDMapping(@NotNull final Map<String, Set<Integer>> doidsPerPrimaryTumor) {
        this.doidsPerPrimaryTumor = doidsPerPrimaryTumor;
    }

    @NotNull
    private static Set<Integer> toIntegerSet(@NotNull String doidsString) {
        if (doidsString.equals(NO_DOIDS_PRESENT)) {
            return Sets.newHashSet();
        }

        Set<Integer> doids = Sets.newHashSet();
        String[] values = doidsString.split(DOID_SEPARATOR);

        for (String value : values) {
            doids.add(Integer.valueOf(value.trim()));
        }
        return doids;
    }

    @Nullable
    public Set<Integer> findDoids(@NotNull String primaryTumorLocation) {
        return doidsPerPrimaryTumor.get(primaryTumorLocation);
    }
}

package com.hartwig.hmftools.common.actionability.cancertype;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class PrimaryTumorToDOIDMapping {

    private static final InputStream TUMOR_LOCATION_MAPPING_RESOURCE =
            PrimaryTumorToDOIDMapping.class.getResourceAsStream("/actionability/primary_tumor_locations_mapping.csv");

    private static final String DOID_SEPARATOR = ";";

    @NotNull
    private final Map<String, Set<String>> doidsPerPrimaryTumor = Maps.newHashMap();

    @NotNull
    public static PrimaryTumorToDOIDMapping createFromResource() throws IOException {
        return new PrimaryTumorToDOIDMapping(TUMOR_LOCATION_MAPPING_RESOURCE);
    }

    private PrimaryTumorToDOIDMapping(@NotNull final InputStream mappingInputStream) throws IOException {
        final CSVParser parser = CSVParser.parse(mappingInputStream, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        for (final CSVRecord record : parser) {
            final String primaryTumorLocation = record.get("primaryTumorLocation");
            final String doids = record.get("doids");

            doidsPerPrimaryTumor.put(primaryTumorLocation, Sets.newHashSet(doids.split(DOID_SEPARATOR)));
        }
    }

    @Nullable
    public Set<String> doidsForPrimaryTumorLocation(@NotNull String primaryTumorLocation) {
        return doidsPerPrimaryTumor.get(primaryTumorLocation);
    }
}

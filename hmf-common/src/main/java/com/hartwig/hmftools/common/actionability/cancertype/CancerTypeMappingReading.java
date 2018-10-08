package com.hartwig.hmftools.common.actionability.cancertype;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CancerTypeMappingReading {

    private static final InputStream TUMOR_LOCATION_MAPPING_RESOURCE =
            CancerTypeMappingReading.class.getResourceAsStream("/actionability/primary_tumor_locations_mapping.csv");

    @NotNull
    private final Map<String, CancerTypeMapping> tumorLocationMap = Maps.newHashMap();

    @NotNull
    public static CancerTypeMappingReading readingFile() throws IOException {
        return new CancerTypeMappingReading(TUMOR_LOCATION_MAPPING_RESOURCE);
    }

    private CancerTypeMappingReading(@NotNull final InputStream mappingInputStream) throws IOException {
        final CSVParser parser = CSVParser.parse(mappingInputStream, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        for (final CSVRecord record : parser) {
            final String primaryTumorLocation = record.get("primaryTumorLocation");
            final String doids = record.get("doids");

            tumorLocationMap.put(primaryTumorLocation,
                    ImmutableCancerTypeMapping.of(primaryTumorLocation, doids));
        }
    }

    @Nullable
    public String doidsForPrimaryTumorLocation(@NotNull String primaryTumorLocation) {
        CancerTypeMapping mapping = tumorLocationMap.get(primaryTumorLocation);
        return mapping != null ? mapping.doids() : null;
    }
}

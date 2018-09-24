package com.hartwig.hmftools.actionability.cancerTypeMapping;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.security.Key;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.jetbrains.annotations.NotNull;

public class CancerTypeMappingReading {
    private static final org.apache.logging.log4j.Logger LOGGER = LogManager.getLogger(CancerTypeMappingReading.class);


    private static final InputStream TUMOR_LOCATION_MAPPING_RESOURCE =
            CancerTypeMappingReading.class.getResourceAsStream("/primary_tumor_locations_mapping.csv");

    @NotNull
    private final Map<String, CancerTypeMapping> tumorLocationMap = Maps.newHashMap();

    @NotNull
    public static CancerTypeMappingReading readingFile() throws IOException {
        return new CancerTypeMappingReading(TUMOR_LOCATION_MAPPING_RESOURCE);
    }

    public CancerTypeMappingReading(@NotNull final InputStream mappingInputStream) throws IOException {
        final CSVParser parser = CSVParser.parse(mappingInputStream, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        for (final CSVRecord record : parser) {
            final String primaryTumorLocation = record.get("primaryTumorLocation");
            final String doids = record.get("doids");

            tumorLocationMap.put(null,
                    ImmutableCancerTypeMapping.of(primaryTumorLocation, doids));
        }
    }
}

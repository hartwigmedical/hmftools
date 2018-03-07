package com.hartwig.hmftools.patientdb.curators;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Map;

import com.google.common.collect.Maps;
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

public class TumorLocationCurator {

    private static final Logger LOGGER = LogManager.getLogger(TumorLocationCurator.class);

    @NotNull
    private final Map<String, CuratedCancerType> tumorLocationMap = Maps.newHashMap();

    public TumorLocationCurator(@NotNull final InputStream mappingInputStream) throws IOException {
        final CSVParser parser = CSVParser.parse(mappingInputStream, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        for (final CSVRecord record : parser) {
            final String location = record.get("primaryTumorLocation");
            final String category = record.get("category");
            final String subcategory = record.get("subcategory");
            tumorLocationMap.put(location.toLowerCase(),
                    ImmutableCuratedCancerType.of(Utils.capitalize(category), Utils.capitalize(subcategory), location));
        }
    }

    @NotNull
    public CuratedCancerType search(@Nullable final String searchTerm) {
        if (searchTerm != null) {
            final CuratedCancerType result = tumorLocationMap.get(searchTerm.toLowerCase());
            if (result != null) {
                return result;
            }
        }

        LOGGER.warn("Could not curate tumor location (using " + System.getProperty("file.encoding") + "): " + searchTerm);
        return ImmutableCuratedCancerType.of(null, null, searchTerm);
    }
}

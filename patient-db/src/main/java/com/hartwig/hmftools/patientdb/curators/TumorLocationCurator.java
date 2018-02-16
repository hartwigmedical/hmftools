package com.hartwig.hmftools.patientdb.curators;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocation;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TumorLocationCurator {
    private final Map<String, CuratedTumorLocation> tumorLocationMap = Maps.newHashMap();

    public TumorLocationCurator(@NotNull final InputStream mappingInputStream) throws IOException {
        final CSVParser parser = CSVParser.parse(mappingInputStream, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        for (final CSVRecord record : parser) {
            final String location = record.get("primaryTumorLocation");
            final String category = record.get("category");
            final String subcategory = record.get("subcategory");
            tumorLocationMap.put(location.toLowerCase(),
                    ImmutableCuratedTumorLocation.of(Utils.capitalize(category), Utils.capitalize(subcategory), location));
        }
    }

    @NotNull
    public CuratedTumorLocation search(@Nullable final String searchTerm) {
        if (searchTerm != null) {
            final CuratedTumorLocation result = tumorLocationMap.get(searchTerm.toLowerCase());
            if (result != null) {
                return result;
            }
        }
        return ImmutableCuratedTumorLocation.of(null, null, searchTerm);
    }
}

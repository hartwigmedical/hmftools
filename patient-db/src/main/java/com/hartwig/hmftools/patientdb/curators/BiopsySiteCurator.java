package com.hartwig.hmftools.patientdb.curators;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.util.Map;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.LoadClinicalData;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.CuratedBiopsyType;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedBiopsyType;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BiopsySiteCurator {

    private static final Logger LOGGER = LogManager.getLogger(BiopsySiteCurator.class);

    private static final InputStream BIOPSY_SITE_MAPPING_RESOURCE = LoadClinicalData.class.getResourceAsStream("/biopsy_site_mapping.csv");

    @NotNull
    private final Map<Key, CuratedBiopsyType> curationMap = Maps.newHashMap();

    @NotNull
    public static BiopsySiteCurator fromProductionResource() throws IOException {
        return new BiopsySiteCurator(BIOPSY_SITE_MAPPING_RESOURCE);
    }

    @VisibleForTesting
    BiopsySiteCurator(@NotNull final InputStream mappingInputStream) throws IOException {
        final CSVParser parser = CSVParser.parse(mappingInputStream, Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        for (final CSVRecord record : parser) {
            final String cancerType = record.get("cancerType");
            final String cancerSubType = record.get("cancerSubType");
            final String biopsySite = record.get("biopsySite");
            final String biopsyLocation = record.get("biopsyLocation");
            final String biopsyType = record.get("biopsyType");

            Key key = Key.fromInputs(cancerType, cancerSubType, biopsySite, biopsyLocation);

            if (curationMap.get(key) != null) {
                LOGGER.warn("Duplicate key found: " + key);
            }

            curationMap.put(key,
                    ImmutableCuratedBiopsyType.of(Utils.capitalize(biopsyType), cancerType, cancerSubType, biopsySite, biopsyLocation));
        }
    }

    @NotNull
    public CuratedBiopsyType search(@Nullable String searchCancerType, @Nullable String searchCancerSubType,
            @Nullable String searchBiopsySite, @Nullable String searchBiopsyLocation) {
        CuratedBiopsyType result = null;
        if (searchCancerType != null && searchCancerSubType != null) {
            result = curationMap.get(Key.fromInputs(searchCancerType, searchCancerSubType, searchBiopsySite, searchBiopsyLocation));
        }

        return result != null
                ? result
                : ImmutableCuratedBiopsyType.of(null, searchCancerType, searchCancerSubType, searchBiopsySite, searchBiopsyLocation);
    }

    private static class Key {
        @NotNull
        private final String cancerType;
        @NotNull
        private final String cancerSubType;
        @Nullable
        private final String biopsySite;
        @Nullable
        private final String biopsyLocation;

        private static Key fromInputs(@NotNull String cancerType, @NotNull String cancerSubType, @Nullable String biopsySite,
                @Nullable String biopsyLocation) {
            String effectiveBiopsySite = biopsySite != null ? biopsySite.toLowerCase() : "null";
            String effectiveBiopsyLocation = biopsyLocation != null ? biopsyLocation.toLowerCase() : "null";

            return new Key(cancerType.toLowerCase(),
                    cancerSubType.toLowerCase(),
                    effectiveBiopsySite.equalsIgnoreCase("null") ? null : effectiveBiopsySite,
                    effectiveBiopsyLocation.equalsIgnoreCase("null") ? null : effectiveBiopsyLocation);
        }

        private Key(@NotNull final String cancerType, @NotNull final String cancerSubType, @Nullable final String biopsySite,
                @Nullable final String biopsyLocation) {
            this.cancerType = cancerType;
            this.cancerSubType = cancerSubType;
            this.biopsySite = biopsySite;
            this.biopsyLocation = biopsyLocation;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final Key key = (Key) o;
            return Objects.equals(cancerType, key.cancerType) && Objects.equals(cancerSubType, key.cancerSubType) && Objects.equals(
                    biopsySite, key.biopsySite) && Objects.equals(biopsyLocation, key.biopsyLocation);
        }

        @Override
        public int hashCode() {
            return Objects.hash(cancerType, cancerSubType, biopsySite, biopsyLocation);
        }

        @Override
        public String toString() {
            return "Key{" + "cancerType='" + cancerType + '\'' + ", cancerSubType='" + cancerSubType + '\'' + ", biopsySite='" + biopsySite
                    + '\'' + ", biopsyLocation='" + biopsyLocation + '\'' + '}';
        }
    }
}

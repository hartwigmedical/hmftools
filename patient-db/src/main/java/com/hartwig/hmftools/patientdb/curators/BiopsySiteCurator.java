package com.hartwig.hmftools.patientdb.curators;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
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

    @NotNull
    private final Map<Key, CuratedBiopsyType> curationMap = Maps.newHashMap();

    @VisibleForTesting
    public BiopsySiteCurator(@NotNull String mappingInput) throws IOException {
        BufferedReader file = new BufferedReader(new FileReader(mappingInput));
        CSVParser parser = CSVFormat.DEFAULT.withDelimiter(',').withHeader().parse(file);

        for (CSVRecord record : parser) {
            String primaryTumorLocation = record.get("primaryTumorLocation");
            String cancerSubType = record.get("cancerSubType");
            String biopsySite = record.get("biopsySite");
            String biopsyLocation = record.get("biopsyLocation");
            String biopsyType = record.get("biopsyType");

            Key key = Key.fromInputs(primaryTumorLocation, cancerSubType, biopsySite, biopsyLocation);

            if (curationMap.get(key) != null) {
                LOGGER.warn("Duplicate key found: '{}'", key);
            }

            curationMap.put(key,
                    ImmutableCuratedBiopsyType.of(Utils.capitalize(biopsyType),
                            primaryTumorLocation,
                            cancerSubType,
                            biopsySite,
                            biopsyLocation));
        }
    }

    @NotNull
    public CuratedBiopsyType search(@Nullable String searchPrimaryTumorLocation, @Nullable String searchCancerSubType,
            @Nullable String searchBiopsySite, @Nullable String searchBiopsyLocation) {
        CuratedBiopsyType result = null;
        if (searchPrimaryTumorLocation != null && searchCancerSubType != null) {
            result = curationMap.get(Key.fromInputs(searchPrimaryTumorLocation,
                    searchCancerSubType,
                    searchBiopsySite,
                    searchBiopsyLocation));
        }

        return result != null
                ? result
                : ImmutableCuratedBiopsyType.of(null,
                        searchPrimaryTumorLocation,
                        searchCancerSubType,
                        searchBiopsySite,
                        searchBiopsyLocation);
    }

    private static class Key {
        @NotNull
        private final String primaryTumorLocation;
        @NotNull
        private final String cancerSubType;
        @Nullable
        private final String biopsySite;
        @Nullable
        private final String biopsyLocation;

        private static Key fromInputs(@NotNull String primaryTumorLocation, @NotNull String cancerSubType, @Nullable String biopsySite,
                @Nullable String biopsyLocation) {
            String effectiveBiopsySite = biopsySite != null ? biopsySite.toLowerCase() : "null";
            String effectiveBiopsyLocation = biopsyLocation != null ? biopsyLocation.toLowerCase() : "null";

            return new Key(primaryTumorLocation.toLowerCase(),
                    cancerSubType.toLowerCase(),
                    effectiveBiopsySite.equalsIgnoreCase("null") ? null : effectiveBiopsySite,
                    effectiveBiopsyLocation.equalsIgnoreCase("null") ? null : effectiveBiopsyLocation);
        }

        private Key(@NotNull final String primaryTumorLocation, @NotNull final String cancerSubType, @Nullable final String biopsySite,
                @Nullable final String biopsyLocation) {
            this.primaryTumorLocation = primaryTumorLocation;
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
            return Objects.equals(primaryTumorLocation, key.primaryTumorLocation) && Objects.equals(cancerSubType, key.cancerSubType)
                    && Objects.equals(biopsySite, key.biopsySite) && Objects.equals(biopsyLocation, key.biopsyLocation);
        }

        @Override
        public int hashCode() {
            return Objects.hash(primaryTumorLocation, cancerSubType, biopsySite, biopsyLocation);
        }

        @Override
        public String toString() {
            return "Key{" + "primaryTumorLocation='" + primaryTumorLocation + '\'' + ", cancerSubType='" + cancerSubType + '\''
                    + ", biopsySite='" + biopsySite + '\'' + ", biopsyLocation='" + biopsyLocation + '\'' + '}';
        }
    }
}

package com.hartwig.hmftools.patientdb.clinical.curators;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedBiopsyType;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableCuratedBiopsyType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BiopsySiteCurator {

    private static final Logger LOGGER = LogManager.getLogger(BiopsySiteCurator.class);
    private static final String FIELD_DELIMITER = "\t";

    @NotNull
    private final Map<Key, CuratedBiopsyType> curationMap = Maps.newHashMap();

    @VisibleForTesting
    public BiopsySiteCurator(@NotNull String biopsyMappingTsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(biopsyMappingTsv).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_DELIMITER);
            String primaryTumorLocation = parts[0];
            String cancerSubType = parts[1];
            String biopsySite = parts[2];
            String biopsyLocation = parts[3];
            String biopsyType = parts[4];

            Key key = Key.fromInputs(primaryTumorLocation, cancerSubType, biopsySite, biopsyLocation);

            if (curationMap.get(key) != null) {
                LOGGER.warn("Duplicate key found: '{}'", key);
            }

            curationMap.put(key,
                    ImmutableCuratedBiopsyType.of(capitalize(biopsyType),
                            primaryTumorLocation,
                            cancerSubType,
                            biopsySite,
                            biopsyLocation));
        }
    }

    @NotNull
    private static String capitalize(@NotNull String string) {
        if (string.isEmpty()) {
            return string;
        } else {
            return string.toUpperCase().substring(0, 1) + string.substring(1);
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

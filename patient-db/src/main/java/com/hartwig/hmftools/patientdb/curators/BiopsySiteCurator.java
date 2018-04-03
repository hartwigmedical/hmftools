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
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class BiopsySiteCurator {

    private static final InputStream BIOPSY_SITE_MAPPING_RESOURCE = LoadClinicalData.class.getResourceAsStream("/biopsy_site_mapping.csv");

    @NotNull
    private final Map<CurationKey, CuratedBiopsyType> biopsyCurationMap = Maps.newHashMap();

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
            biopsyCurationMap.put(new CurationKey(cancerType.toLowerCase(),
                            cancerSubType.toLowerCase(),
                            biopsySite.toLowerCase(),
                            biopsyLocation.toLowerCase()),
                    ImmutableCuratedBiopsyType.of(Utils.capitalize(biopsyType), cancerType, cancerSubType, biopsySite, biopsyLocation));
        }
    }

    @NotNull
    public CuratedBiopsyType search(@Nullable String searchCancerType, @Nullable String searchCancerSubType,
            @Nullable String searchBiopsySite, @Nullable String searchBiopsyLocation) {
        if (searchCancerType != null && searchCancerSubType != null && searchBiopsySite != null && searchBiopsyLocation != null) {
            CurationKey effectiveSearchTerm = new CurationKey(searchCancerType.toLowerCase(),
                    searchCancerSubType.toLowerCase(),
                    searchBiopsySite.toLowerCase(),
                    searchBiopsyLocation.toLowerCase());
            final CuratedBiopsyType result = biopsyCurationMap.get(effectiveSearchTerm);

            if (result != null) {
                return result;
            }
        }

        return ImmutableCuratedBiopsyType.of(null, searchCancerType, searchCancerSubType, searchBiopsySite, searchBiopsyLocation);
    }

    private static class CurationKey {
        @NotNull
        private final String cancerType;
        @NotNull
        private final String cancerSubType;
        @NotNull
        private final String biopsySite;
        @NotNull
        private final String biopsyLocation;

        private CurationKey(@NotNull final String cancerType, @NotNull final String cancerSubType, @NotNull final String biopsySite,
                @NotNull final String biopsyLocation) {
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
            final CurationKey that = (CurationKey) o;
            return Objects.equals(cancerType, that.cancerType) && Objects.equals(cancerSubType, that.cancerSubType) && Objects.equals(
                    biopsySite,
                    that.biopsySite) && Objects.equals(biopsyLocation, that.biopsyLocation);
        }

        @Override
        public int hashCode() {

            return Objects.hash(cancerType, cancerSubType, biopsySite, biopsyLocation);
        }
    }
}

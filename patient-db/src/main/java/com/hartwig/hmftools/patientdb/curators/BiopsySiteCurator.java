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
    private final Map<CancerTypeBiopsySite, CuratedBiopsyType> biopsySiteMap = Maps.newHashMap();

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
            final String type = record.get("type");
            biopsySiteMap.put(new CancerTypeBiopsySite(cancerType.toLowerCase(), cancerSubType.toLowerCase(), biopsySite.toLowerCase()),
                    ImmutableCuratedBiopsyType.of(Utils.capitalize(type), cancerType, cancerSubType, biopsySite));
        }
    }

    @NotNull
    public CuratedBiopsyType search(@Nullable String searchCancerType, @Nullable String searchCancerSubType,
            @Nullable String searchBiopsySite) {
        if (searchCancerType != null && searchBiopsySite != null && searchCancerSubType != null) {
            CancerTypeBiopsySite effectiveSearchTerm = new CancerTypeBiopsySite(searchCancerType.toLowerCase(),
                    searchCancerSubType.toLowerCase(),
                    searchBiopsySite.toLowerCase());
            final CuratedBiopsyType result = biopsySiteMap.get(effectiveSearchTerm);

            if (result != null) {
                return result;
            }
        }

        return ImmutableCuratedBiopsyType.of(null, searchCancerType, searchCancerSubType, searchBiopsySite);
    }

    private static class CancerTypeBiopsySite {
        @NotNull
        private final String cancerType;
        @NotNull
        private final String cancerSubType;
        @NotNull
        private final String biopsySite;

        private CancerTypeBiopsySite(@NotNull final String cancerType, @NotNull final String cancerSubType,
                @NotNull final String biopsySite) {
            this.cancerType = cancerType;
            this.cancerSubType = cancerSubType;
            this.biopsySite = biopsySite;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final CancerTypeBiopsySite that = (CancerTypeBiopsySite) o;
            return Objects.equals(cancerType, that.cancerType) && Objects.equals(cancerSubType, that.cancerSubType) && Objects.equals(
                    biopsySite,
                    that.biopsySite);
        }

        @Override
        public int hashCode() {
            return Objects.hash(cancerType, cancerSubType, biopsySite);
        }
    }
}

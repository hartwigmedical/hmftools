package com.hartwig.hmftools.patientdb.curators;

import java.io.IOException;

import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;

public final class TestCuratorFactory {

    public static final String TUMOR_LOCATION_MAPPING_CSV = Resources.getResource("curators/test_tumor_location_mapping.csv").getPath();
    public static final String TUMOR_LOCATION_V2_MAPPING_TSV =
            Resources.getResource("curators/test_tumor_location_v2_mapping.tsv").getPath();
    public static final String BIOPSY_SITE_MAPPING_CSV = Resources.getResource("curators/test_biopsy_site_mapping.csv").getPath();
    public static final String TREATMENT_MAPPING_CSV = Resources.getResource("curators/test_treatment_mapping.csv").getPath();

    private TestCuratorFactory() {
    }

    @NotNull
    public static TumorLocationCuratorV2 tumorLocationV2Curator() {
        try {
            return new TumorLocationCuratorV2(TUMOR_LOCATION_V2_MAPPING_TSV);
        } catch (IOException e) {
            throw new IllegalStateException("Could not create tumor location V2 curator!");
        }
    }

    @NotNull
    public static TumorLocationCurator tumorLocationCurator() {
        try {
            return new TumorLocationCurator(TUMOR_LOCATION_MAPPING_CSV);
        } catch (IOException e) {
            throw new IllegalStateException("Could not create tumor location curator!");
        }
    }

    @NotNull
    public static BiopsySiteCurator biopsySiteCurator() {
        try {
            return new BiopsySiteCurator(BIOPSY_SITE_MAPPING_CSV);
        } catch (IOException e) {
            throw new IllegalStateException("Could not create biopsy site curator!");
        }
    }

    @NotNull
    public static TreatmentCurator treatmentCurator() {
        try {
            return new TreatmentCurator(TREATMENT_MAPPING_CSV);
        } catch (IOException e) {
            throw new IllegalStateException("Could not create treatment curator!");
        }
    }
}

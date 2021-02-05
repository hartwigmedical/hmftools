package com.hartwig.hmftools.patientdb.clinical;

public final class ClinicalConstants {

    private ClinicalConstants() {
    }

    public static final int MAX_DAYS_BETWEEN_SAMPLING_AND_BIOPSY_DATE = 30;
    public static final int MAX_DAYS_ARRIVAL_DATE_AFTER_BIOPSY_DATE = 180;
    public static final int MAX_DAYS_BETWEEN_TREATMENT_AND_BIOPSY = 90;
    public static final int START_DATE_RESPONSE_THRESHOLD = 112;
    public static final int IMMEDIATE_TREATMENT_END_THRESHOLD = 7;

}

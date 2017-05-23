package com.hartwig.hmftools.patientdb;

public final class Config {

    private Config() {
    }

    public static final int SAMPLING_DATE_THRESHOLD = 30;
    public static final int ARRIVAL_DATE_THRESHOLD = 180;
    public static final int MAX_DAYS_BETWEEN_TREATMENT_AND_BIOPSY = 90;
}

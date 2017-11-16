package com.hartwig.hmftools.patientdb;

public final class Config {

    private Config() {
    }

    public static final int MAX_BIOPSY_DAYS_PRIOR_TO_REG_DATE = 30;
    public static final int MAX_DAYS_BETWEEN_SAMPLING_AND_BIOPSY_DATE = 30;
    public static final int ARRIVAL_DATE_THRESHOLD = 180;
    public static final int MAX_DAYS_BETWEEN_TREATMENT_AND_BIOPSY = 90;
    public static final int BATCH_INSERT_SIZE = 1000;
    public static final int START_DATE_RESPONSE_THRESHOLD = 112;
    public static final int IMMEDIATE_TREATMENT_END_THRESHOLD = 7;
}

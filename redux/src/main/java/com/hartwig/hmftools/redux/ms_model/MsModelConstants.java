package com.hartwig.hmftools.redux.ms_model;

import static com.hartwig.hmftools.redux.ms_sites.UnitKey.A_T;

import java.util.List;

public final class MsModelConstants
{
    public static final String REPEAT_3_5_GROUP = "3-5bp repeat";
    public static final String MULTI_BASE_REPEAT = "bp repeat";

    // training constants
    public static final int TRAINING_NOISE_THRESHOLD_DEFAULT = 90;
    public static final double TRAINING_SAMPLE_PURITY_THRESHOLD_DEFAULT = 0.5;

    public static final int TRAINING_REPEAT_COUNT_MIN_DEFAULT = 10;
    public static final int TRAINING_REPEAT_COUNT_MAX_DEFAULT = 20;
    public static final List<String> TRAINING_REPEAT_UNITS_DEFAULT = List.of(A_T.keyStr());
}

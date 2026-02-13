package com.hartwig.hmftools.redux.ms_model;

import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_NOISE_THRESHOLD_DEFAULT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_REPEAT_COUNT_MAX_DEFAULT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_REPEAT_COUNT_MIN_DEFAULT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_REPEAT_UNITS_DEFAULT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_SAMPLE_PURITY_THRESHOLD_DEFAULT;

import java.util.List;

public class MsModelParams
{
    public final int NoiseThreshold;
    public final double SamplePurityThreshold;

    public final List<String> RepeatTypes;
    public final int RepeatCountMin;
    public final int RepeatCountMax;

    public MsModelParams(
            final int noiseThreshold, final double samplePurityThreshold, final List<String> repeatTypes,
            final int repeatCountMin, final int repeatCountMax)
    {
        NoiseThreshold = noiseThreshold;
        SamplePurityThreshold = samplePurityThreshold;
        RepeatTypes = repeatTypes;
        RepeatCountMin = repeatCountMin;
        RepeatCountMax = repeatCountMax;
    }

    public static final MsModelParams DEFAULT_MODEL_PARAMS = new MsModelParams(
            TRAINING_NOISE_THRESHOLD_DEFAULT,
            TRAINING_SAMPLE_PURITY_THRESHOLD_DEFAULT,
            TRAINING_REPEAT_UNITS_DEFAULT,
            TRAINING_REPEAT_COUNT_MIN_DEFAULT,
            TRAINING_REPEAT_COUNT_MAX_DEFAULT);

}

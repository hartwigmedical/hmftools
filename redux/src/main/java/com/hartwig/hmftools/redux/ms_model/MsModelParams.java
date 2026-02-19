package com.hartwig.hmftools.redux.ms_model;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_NOISE_THRESHOLD_DEFAULT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_REPEAT_COUNT_MAX_DEFAULT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_REPEAT_COUNT_MIN_DEFAULT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_REPEAT_UNITS_DEFAULT;
import static com.hartwig.hmftools.redux.ms_model.MsModelConstants.TRAINING_SAMPLE_PURITY_THRESHOLD_DEFAULT;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class MsModelParams
{
    public final int NoiseThreshold;
    public final double SamplePurityThreshold;

    public final List<String> RepeatTypes;
    public final int RepeatCountMin;
    public final int RepeatCountMax;

    public static final String ERROR_RATE_NOISE_THRESHOLD = "error_rate_noise";
    public static final String PURITY_THRESHOLD = "purity_threshold";
    public static final String FITTING_UNITS = "fitting_repeat_units";
    public static final String FITTING_COUNT_MIN = "fitting_count_min";
    public static final String FITTING_COUNT_MAX = "fitting_count_max";

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

    public static MsModelParams fromConfig(final ConfigBuilder configBuilder)
    {
        int noiseThreshold = configBuilder.getInteger(ERROR_RATE_NOISE_THRESHOLD);
        double samplePurityThreshold = configBuilder.getDecimal(PURITY_THRESHOLD);

        List<String> repeatTypes = Lists.newArrayList();

        if(configBuilder.hasValue(FITTING_UNITS))
        {
            String[] repeatUnits = configBuilder.getValue(FITTING_UNITS).split(ITEM_DELIM, -1);
            Arrays.stream(repeatUnits).forEach(x -> repeatTypes.add(x));
        }
        else
        {
            repeatTypes.addAll(TRAINING_REPEAT_UNITS_DEFAULT);
        }

        int repeatCountMin = configBuilder.getInteger(FITTING_COUNT_MIN);
        int repeatCountMax = configBuilder.getInteger(FITTING_COUNT_MAX);

        return new MsModelParams(noiseThreshold, samplePurityThreshold, repeatTypes, repeatCountMin, repeatCountMax);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(ERROR_RATE_NOISE_THRESHOLD, "Error rate noise threshold", TRAINING_NOISE_THRESHOLD_DEFAULT);
        configBuilder.addDecimal(PURITY_THRESHOLD, "Sample purity threshold", TRAINING_SAMPLE_PURITY_THRESHOLD_DEFAULT);
        configBuilder.addConfigItem(FITTING_UNITS, "Repeat units separated by ';'");
        configBuilder.addInteger(FITTING_COUNT_MIN, "Repeat unit count minimum", TRAINING_REPEAT_COUNT_MIN_DEFAULT);
        configBuilder.addInteger(FITTING_COUNT_MAX, "Repeat unit count maximum", TRAINING_REPEAT_COUNT_MAX_DEFAULT);
    }

    public static final MsModelParams DEFAULT_MODEL_PARAMS = new MsModelParams(
            TRAINING_NOISE_THRESHOLD_DEFAULT,
            TRAINING_SAMPLE_PURITY_THRESHOLD_DEFAULT,
            TRAINING_REPEAT_UNITS_DEFAULT,
            TRAINING_REPEAT_COUNT_MIN_DEFAULT,
            TRAINING_REPEAT_COUNT_MAX_DEFAULT);
}

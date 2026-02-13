package com.hartwig.hmftools.redux.ms_model;

import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;

public class MsModelConfig
{
    public final List<String> SampleIds;
    public final String ReduxDir;
    public final String PurpleDir;
    public final String OutputDir;
    public final String OutputId;
    public final String ModelCoefficientsFile;
    public final String ModelErroRatesFile;

    public final List<TrainingRoutines> Routines;

    private static final String ROUTINES = "routines";

    private static final String MODEL_COEFF_FILE = "model_coefficients";
    private static final String MODEL_ERROR_RATES_FILE = "model_error_rates";

    protected enum TrainingRoutines
    {
        COEFFICIENTS,
        ERROR_RATES,
        VALIDATION;
    }

    public MsModelConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = Lists.newArrayList();
        SampleIds.addAll(loadSampleIdsFile(configBuilder));

        Routines = Lists.newArrayList();
        String modesStr = configBuilder.getValue(ROUTINES);

        if(modesStr.equals("ALL"))
        {
            Arrays.stream(TrainingRoutines.values()).forEach(x -> Routines.add(x));
        }
        else
        {
            Arrays.stream(modesStr.split(ITEM_DELIM, -1)).forEach(x -> Routines.add(TrainingRoutines.valueOf(x)));
        }

        ModelCoefficientsFile = configBuilder.getValue(MODEL_COEFF_FILE);
        ModelErroRatesFile = configBuilder.getValue(MODEL_ERROR_RATES_FILE);

        PurpleDir = checkAddDirSeparator(configBuilder.getValue(PURPLE_DIR_CFG));
        ReduxDir = checkAddDirSeparator(configBuilder.getValue(REDUX_DIR_CFG));
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        addSampleIdFile(configBuilder, true);
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addPath(REDUX_DIR_CFG, true, REDUX_DIR_DESC);
        configBuilder.addPath(MODEL_COEFF_FILE, false, "Model coefficients file");
        configBuilder.addPath(MODEL_ERROR_RATES_FILE, false, "Model cohort error rates file");

        configBuilder.addConfigItem(ROUTINES, true, enumValueSelectionAsStr(TrainingRoutines.values(), "Routines"));

        ConfigUtils.addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);
    }
}

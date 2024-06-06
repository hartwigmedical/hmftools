package com.hartwig.hmftools.cup.cli;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class PredictionConfig
{
    public final String SampleId;
    public final String ClassifierPath;
    public final String FeaturesPath;
    public final String OutputDir;
    public final String PythonPath;

    public static final String CLASSIFIER_PATH = "classifier_path";
    public static final String CLASSIFIER_PATH_DESC = "Path to the CUPPA classifier file (.pickle or .pickle.gz file)";
    public static final String FEATURES_PATH = "features_path";
    public static final String FEATURES_PATH_DESC = "Path to or directory containing the input features file(s) (.cuppa_data.tsv file(s))";
    public static final String PYTHON_PATH = "python_path";
    public static final String PYTHON_PATH_DESC = "Path to Python interpreter binary";

    public PredictionConfig(ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        ClassifierPath = configBuilder.getValue(CLASSIFIER_PATH);
        FeaturesPath = configBuilder.getValue(FEATURES_PATH, "");
        OutputDir = parseOutputDir(configBuilder);
        PythonPath = configBuilder.getValue(PYTHON_PATH);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(FEATURES_PATH, false, CLASSIFIER_PATH_DESC);
        configBuilder.addPath(CLASSIFIER_PATH, true, FEATURES_PATH_DESC);
        configBuilder.addPath(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        configBuilder.addPath(PYTHON_PATH, true, PYTHON_PATH_DESC);
    }
}

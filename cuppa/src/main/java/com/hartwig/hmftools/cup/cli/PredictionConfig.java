package com.hartwig.hmftools.cup.cli;

import com.hartwig.hmftools.common.utils.config.CommonConfig;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

public class PredictionConfig
{
    public final String SampleId;
    public final String ClassifierPath;
    public final String FeaturesPath;
    public final String OutputDir;
    public final String PythonPath;

    public final String ClfGroup;
    public final String CvPredictionsPath;

    public static final String CLASSIFIER_PATH = "classifier_path";
    public static final String CLASSIFIER_PATH_DESC = "Path to the CUPPA classifier file (.pickle or .pickle.gz file)";

    public static final String FEATURES_PATH = "features_path";
    public static final String FEATURES_PATH_DESC = "Path to or directory containing the input features file(s), i.e. *.cuppa_data.tsv file(s)";

    public static final String CLF_GROUP = "clf_group";
    public static final String CLF_GROUP_DESC = "Classifier groups to filter probabilities for. Can be 'all' (default) or 'dna'";

    public static final String CV_PREDICTIONS_PATH = "cv_predictions_path";
    public static final String CV_PREDICTIONS_PATH_DESC = "Path to a CuppaPrediction tsv file containing the cross-validation predictions";

    public static final String PYTHON_PATH = "python_path";
    public static final String PYTHON_PATH_DESC = "Path to Python interpreter binary";

    public PredictionConfig(ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(CommonConfig.SAMPLE);
        ClassifierPath = configBuilder.getValue(CLASSIFIER_PATH);
        FeaturesPath = configBuilder.getValue(FEATURES_PATH);
        OutputDir = FileWriterUtils.parseOutputDir(configBuilder);
        PythonPath = configBuilder.getValue(PYTHON_PATH);

        ClfGroup = configBuilder.getValue(CLF_GROUP);
        CvPredictionsPath = configBuilder.getValue(CV_PREDICTIONS_PATH);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(CommonConfig.SAMPLE, false, CommonConfig.SAMPLE_DESC);
        configBuilder.addPath(FEATURES_PATH, false, FEATURES_PATH_DESC);
        configBuilder.addPath(CLASSIFIER_PATH, true, CLASSIFIER_PATH_DESC);
        configBuilder.addPath(FileWriterUtils.OUTPUT_DIR, true, FileWriterUtils.OUTPUT_DIR_DESC);
        configBuilder.addPath(PYTHON_PATH, true, PYTHON_PATH_DESC);

        configBuilder.addConfigItem(CLF_GROUP, false, CLF_GROUP_DESC);
        configBuilder.addConfigItem(CV_PREDICTIONS_PATH, false, CV_PREDICTIONS_PATH_DESC);
    }
}

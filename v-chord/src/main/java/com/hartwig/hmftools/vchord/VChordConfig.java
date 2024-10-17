package com.hartwig.hmftools.vchord;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class VChordConfig
{
    private static final String MODEL_PATH_CFG = "model";

    private final String mSampleId;
    private final String mPurpleDir;
    private final String mModelPath;

    private final String mOutputDir;

    public String getSampleId()
    {
        return mSampleId;
    }

    public String getPurpleDir()
    {
        return mPurpleDir;
    }

    public String getModelPath()
    {
        return mModelPath;
    }

    public String getOutputDir()
    {
        return mOutputDir;
    }

    public VChordConfig(final ConfigBuilder configBuilder)
    {
        mSampleId = configBuilder.getValue(SAMPLE);
        mPurpleDir = checkAddDirSeparator(configBuilder.getValue(PURPLE_DIR_CFG));
        mModelPath = configBuilder.getValue(MODEL_PATH_CFG);
        mOutputDir = parseOutputDir(configBuilder);
        checkCreateOutputDir(mOutputDir);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addPath(MODEL_PATH_CFG, true, "path to the torchscript model");
        addOutputDir(configBuilder, true);
        addLoggingOptions(configBuilder);
    }
}
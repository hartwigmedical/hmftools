package com.hartwig.hmftools.orange;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.orange.util.PathResolver;

public class OrangeRnaConfig
{
    public final String RnaSampleId;
    public final String IsofoxDir;

    private static String RNA_SAMPLE_ID = "rna_sample_id";

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(RNA_SAMPLE_ID, false, "(Optional) The RNA sample of the tumor sample for which ORANGE will run");
        configBuilder.addPath(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);
    }

    public OrangeRnaConfig(
            final ConfigBuilder configBuilder, final String tumorSampleId, final PathResolver pathResolver,
            final PipelineToolDirectories defaultToolDirectories)
    {
        if(!configBuilder.hasValue(RNA_SAMPLE_ID) || !configBuilder.hasValue(ISOFOX_DIR_CFG))
        {
            RnaSampleId = null;
            IsofoxDir = null;

            LOGGER.info("RNA config not present, will continue without RNA configuration");
            return;
        }

        RnaSampleId = configBuilder.getValue(RNA_SAMPLE_ID);
        IsofoxDir = pathResolver.resolveMandatoryToolDirectory(ISOFOX_DIR_CFG, defaultToolDirectories.isofoxDir());

        LOGGER.debug("RNA sample configured as {}", RnaSampleId);
    }

    public OrangeRnaConfig(final String rnaSampleId, final String isofoxDir)
    {
        RnaSampleId = rnaSampleId;
        IsofoxDir = isofoxDir;
    }
}

package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.SAMPLE_PROBES_MAX_DEFAULT;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public record SampleVariantsConfig(
        String sampleId,
        String purpleDir,
        @Nullable String linxDir,
        @Nullable String linxGermlineDir,
        int maxProbes
)
{
    private static final String CFG_MAX_PROBES = "sample_probes";
    private static final String DESC_MAX_PROBES = "Maximum number of sample variant probes";

    private static final Logger LOGGER = LogManager.getLogger(SampleVariantsConfig.class);

    @Nullable
    public static SampleVariantsConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        String sampleId = configBuilder.getValue(SAMPLE);
        String purpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        String linxDir = configBuilder.getValue(LINX_DIR_CFG);
        String linxGermlineDir = configBuilder.getValue(LINX_GERMLINE_DIR_CFG);
        int maxProbes = configBuilder.getInteger(CFG_MAX_PROBES);

        // If sampleId is present then assume the user wants to generate sample variant probes, otherwise no.
        if(sampleId == null)
        {
            return null;
        }
        else
        {
            if(purpleDir == null)
            {
                LOGGER.error("Required: {}", PURPLE_DIR_CFG);
                System.exit(1);
            }
            return new SampleVariantsConfig(sampleId, purpleDir, linxDir, linxGermlineDir, maxProbes);
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_DESC);
        configBuilder.addInteger(CFG_MAX_PROBES, DESC_MAX_PROBES, SAMPLE_PROBES_MAX_DEFAULT);
    }
}

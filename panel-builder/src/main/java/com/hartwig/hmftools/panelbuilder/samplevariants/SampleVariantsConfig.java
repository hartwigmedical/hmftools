package com.hartwig.hmftools.panelbuilder.samplevariants;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

public record SampleVariantsConfig(
        String sampleId,
        String purpleDir,
        @Nullable String linxDir,
        @Nullable String linxGermlineDir,
        @Nullable String referenceVariantsFile
)
{
    private static final String CFG_REF_VARIANTS_FILE = "ref_variants";
    private static final String DESC_REF_VARIANTS_FILE = "Reference variants file";

    private static final Logger LOGGER = LogManager.getLogger(SampleVariantsConfig.class);

    @Nullable
    public static SampleVariantsConfig fromConfigBuilder(final ConfigBuilder configBuilder)
    {
        String sampleId = configBuilder.getValue(SAMPLE);
        String purpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        String linxDir = configBuilder.getValue(LINX_DIR_CFG);
        String linxGermlineDir = configBuilder.getValue(LINX_GERMLINE_DIR_CFG);
        String refVariants = configBuilder.getValue(CFG_REF_VARIANTS_FILE);

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
            return new SampleVariantsConfig(sampleId, purpleDir, linxDir, linxGermlineDir, refVariants);
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_DESC);
        configBuilder.addPath(CFG_REF_VARIANTS_FILE, false, DESC_REF_VARIANTS_FILE);
    }
}

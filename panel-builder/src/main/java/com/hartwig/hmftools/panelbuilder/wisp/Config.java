package com.hartwig.hmftools.panelbuilder.wisp;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.Nullable;

public record Config(
        String sampleId,
        String purpleDir,
        @Nullable
        String linxDir,
        @Nullable
        String linxGermlineDir,
        @Nullable
        String referenceVariantsFile
)
{
    private static final String CFG_REF_VARIANTS_FILE = "ref_variants";
    private static final String DESC_REF_VARIANTS_FILE = "Reference variants file";

    public Config(final ConfigBuilder configBuilder)
    {
        this(
                configBuilder.getValue(SAMPLE),
                configBuilder.getValue(PURPLE_DIR_CFG),
                configBuilder.getValue(LINX_DIR_CFG),
                configBuilder.getValue(LINX_GERMLINE_DIR_CFG),
                configBuilder.getValue(CFG_REF_VARIANTS_FILE));
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(LINX_GERMLINE_DIR_CFG, false, LINX_GERMLINE_DIR_DESC);
        configBuilder.addPath(CFG_REF_VARIANTS_FILE, false, DESC_REF_VARIANTS_FILE);
    }
}

package com.hartwig.hmftools.sage.append;

import static com.hartwig.hmftools.sage.quality.QualityConfig.HIGH_DEPTH_MODE;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.sage.SageConfig;

public class SageAppendConfig
{
    public final SageConfig Common;
    public final String InputVcf;
    public final boolean FilterToGenes;

    private static final String INPUT_VCF = "input_vcf";
    private static final String FILTER_TO_GENES = "require_gene";

    public SageAppendConfig(final String version, final ConfigBuilder configBuilder)
    {
        Common = new SageConfig(version, configBuilder, true);

        InputVcf = configBuilder.getValue(INPUT_VCF);
        FilterToGenes = configBuilder.hasFlag(FILTER_TO_GENES);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        SageConfig.registerCommonConfig(configBuilder);
        configBuilder.addPath(INPUT_VCF, true, "Path to input vcf");
        configBuilder.addFlag(FILTER_TO_GENES, "Only process variants with gene annotations");
    }
}

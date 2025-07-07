package com.hartwig.hmftools.sage.vis;

import static com.hartwig.hmftools.common.utils.config.ConfigItemType.INTEGER;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.SimpleVariant;

public class VisConfig
{
    public final boolean Enabled;
    public final boolean PassOnly;
    public final int MaxSupportReads;
    public final String OutputDir;
    public final List<SimpleVariant> SpecificVariants;

    private static final String VIS_OUTPUT_DIR = "vis_output_dir";
    private static final String SPECIFIC_VARIANTS = "vis_variants";
    private static final String PASS_ONLY = "vis_pass_only";
    private static final String MAX_SUPPORT_READS = "vis_max_support_reads";

    private static final String DEFAULT_PLOT_DIR = "vis";

    public VisConfig(final ConfigBuilder configBuilder, final String outputDir)
    {
        PassOnly = configBuilder.hasFlag(PASS_ONLY);
        SpecificVariants = Lists.newArrayList();

        if(configBuilder.hasValue(SPECIFIC_VARIANTS))
            SpecificVariants.addAll(SimpleVariant.fromConfig(configBuilder.getValue(SPECIFIC_VARIANTS)));

        boolean enabled = PassOnly || !SpecificVariants.isEmpty();

        String visDir = null;

        if(configBuilder.hasValue(VIS_OUTPUT_DIR))
        {
            enabled = true;
            visDir = checkAddDirSeparator(outputDir + configBuilder.getValue(VIS_OUTPUT_DIR));
        }
        else if(enabled)
        {
            visDir = outputDir + DEFAULT_PLOT_DIR + File.separator;
        }

        OutputDir = visDir;
        Enabled = enabled;

        MaxSupportReads = configBuilder.getInteger(MAX_SUPPORT_READS);
    }

    public boolean processVariant(final SimpleVariant variant)
    {
        if(SpecificVariants.isEmpty())
            return true;

        return SpecificVariants.stream().anyMatch(x -> variant.matches(x));
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(VIS_OUTPUT_DIR, "Visualiser: HTML output dir");

        configBuilder.addFlag(PASS_ONLY, "Visualiser: write files for passing variants");

        configBuilder.addConfigItem(SPECIFIC_VARIANTS,
                "Visualiser: write files for specific variants. Format: chromosome:position:ref:alt separated by ';'");

        configBuilder.addConfigItem(
                INTEGER, MAX_SUPPORT_READS, false,
                "Visualiser: Max reads by support type, default is fixed by type. Use '-1' to show all.", "0");
    }

    public VisConfig()
    {
        OutputDir = null;
        SpecificVariants = Lists.newArrayList();
        PassOnly = false;
        MaxSupportReads = 0;
        Enabled = false;
    }
}

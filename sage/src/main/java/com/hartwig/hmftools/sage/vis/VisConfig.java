package com.hartwig.hmftools.sage.vis;

import static com.hartwig.hmftools.common.utils.config.ConfigItemType.INTEGER;
import static com.hartwig.hmftools.sage.vis.SageVisConstants.REPORTED_KEY;

import java.io.File;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.variantcontext.VariantContext;

public class VisConfig
{
    public final boolean Enabled;
    public final boolean PassOnly;
    public final int MaxSupportReads;
    public final File OutputDir;
    public final List<SimpleVariant> SpecificVariants;
    public final File PurpleVcf;

    public boolean mNoVariants;

    private static final String VIS_OUTPUT_DIR = "vis_output_dir";
    private static final String SPECIFIC_VARIANTS = "vis_variants";
    private static final String PASS_ONLY = "vis_pass_only";
    private static final String MAX_SUPPORT_READS = "vis_max_support_reads";
    private static final String PURPLE_VCF = "vis_purple_vcf";

    private static final String DEFAULT_PLOT_DIR = "vis";

    public VisConfig(final ConfigBuilder configBuilder, final String outputDir)
    {
        PassOnly = configBuilder.hasFlag(PASS_ONLY);
        SpecificVariants = Lists.newArrayList();

        String VcfStr = configBuilder.hasValue(PURPLE_VCF) ? configBuilder.getValue(PURPLE_VCF) : null;
        PurpleVcf = VcfStr == null ? null : new File(VcfStr);

        mNoVariants = false;
        if(configBuilder.hasValue(SPECIFIC_VARIANTS))
        {
            SpecificVariants.addAll(SimpleVariant.fromConfig(configBuilder.getValue(SPECIFIC_VARIANTS)));
        }
        else if(PurpleVcf != null)
        {
            try(VcfFileReader vcfFileReader = new VcfFileReader(PurpleVcf.toString(), true))
            {
                CloseableTribbleIterator<VariantContext> iter = vcfFileReader.iterator();
                while(iter.hasNext())
                {
                    VariantContext variant = iter.next();
                    Object reported = variant.getAttribute(REPORTED_KEY);
                    if(reported != null && (boolean) reported)
                    {
                        String contig = variant.getContig();
                        int pos = variant.getStart();
                        String ref = variant.getReference().getBaseString();
                        String alt = variant.getAltAlleleWithHighestAlleleCount().getBaseString();
                        SpecificVariants.add(new SimpleVariant(contig, pos, ref, alt));
                    }
                }

                iter.close();
            }

            if(SpecificVariants.isEmpty())
                mNoVariants = true;
        }

        boolean enabled = PassOnly || !SpecificVariants.isEmpty();
        File visDir = null;
        if(configBuilder.hasValue(VIS_OUTPUT_DIR))
        {
            enabled = true;
            visDir = new File(outputDir, configBuilder.getValue(VIS_OUTPUT_DIR));
        }
        else if(enabled)
        {
            visDir = new File(outputDir, DEFAULT_PLOT_DIR);
        }

        OutputDir = visDir;
        Enabled = enabled;

        MaxSupportReads = configBuilder.getInteger(MAX_SUPPORT_READS);
    }

    public boolean processVariant(final SimpleVariant variant)
    {
        if(SpecificVariants.isEmpty())
            return true;

        return SpecificVariants.stream().anyMatch(variant::matches);
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

        configBuilder.addConfigItem(PURPLE_VCF, false, "VCF file containing pave and purple annotations");
    }

    public VisConfig()
    {
        OutputDir = null;
        SpecificVariants = Lists.newArrayList();
        PassOnly = false;
        MaxSupportReads = 0;
        Enabled = false;
        mNoVariants = false;
        PurpleVcf = null;
    }

    public boolean noVariants() { return mNoVariants; }
}

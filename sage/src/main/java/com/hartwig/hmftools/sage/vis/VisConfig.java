package com.hartwig.hmftools.sage.vis;

import static com.hartwig.hmftools.common.utils.config.ConfigItemType.INTEGER;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.SimpleVariant;
import com.hartwig.hmftools.common.variant.VcfFileReader;

import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.variant.variantcontext.VariantContext;

public class VisConfig
{
    public final boolean PassOnly;
    public final int MaxSupportReads;
    public final String OutputDir;

    public final List<SimpleVariant> SpecificVariants;

    private final String mPurpleVcf;
    private final List<VariantAndContext> mReportableVariants;

    private static final String VIS_OUTPUT_DIR = "vis_output_dir";
    private static final String SPECIFIC_VARIANTS = "vis_variants";
    private static final String PASS_ONLY = "vis_pass_only";
    private static final String MAX_SUPPORT_READS = "vis_max_support_reads";
    private static final String PURPLE_VCF = "vis_purple_vcf";

    private static final String DEFAULT_PLOT_DIR = "vis";

    private class VariantAndContext
    {
        public final SimpleVariant Variant;
        public final VariantContext Context;

        public VariantAndContext(final SimpleVariant variant, final VariantContext context)
        {
            Variant = variant;
            Context = context;
        }
    }

    public VisConfig(final ConfigBuilder configBuilder, final String outputDir)
    {
        PassOnly = configBuilder.hasFlag(PASS_ONLY);
        SpecificVariants = Lists.newArrayList();
        MaxSupportReads = configBuilder.getInteger(MAX_SUPPORT_READS);

        mReportableVariants = Lists.newArrayList();
        mPurpleVcf = configBuilder.getValue(PURPLE_VCF);

        if(configBuilder.hasValue(SPECIFIC_VARIANTS))
        {
            SpecificVariants.addAll(SimpleVariant.fromConfig(configBuilder.getValue(SPECIFIC_VARIANTS)));

            SG_LOGGER.info("loaded {} specific variants for vis plot generation", SpecificVariants.size());
        }
        else if(mPurpleVcf != null)
        {
            try(VcfFileReader vcfFileReader = new VcfFileReader(mPurpleVcf, true))
            {
                CloseableTribbleIterator<VariantContext> iter = vcfFileReader.iterator();
                while(iter.hasNext())
                {
                    VariantContext variantContext = iter.next();

                    if(variantContext.getAttributeAsBoolean(REPORTED_FLAG, false))
                    {
                        String chromosome = variantContext.getContig();
                        int pos = variantContext.getStart();
                        String ref = variantContext.getReference().getBaseString();
                        String alt = variantContext.getAltAlleleWithHighestAlleleCount().getBaseString();

                        SimpleVariant variant = new SimpleVariant(chromosome, pos, ref, alt);
                        SpecificVariants.add(variant);
                        mReportableVariants.add(new VariantAndContext(variant, variantContext));
                    }
                }

                iter.close();
            }

            SG_LOGGER.info("loaded {} reportable variants for vis plot generation", SpecificVariants.size());
        }

        String defaultPlotDir = checkAddDirSeparator(outputDir) + DEFAULT_PLOT_DIR;
        OutputDir = configBuilder.getValue(VIS_OUTPUT_DIR, defaultPlotDir);

        // always make the output directory if no vis variants will be written
        if(enabled())
            checkCreateOutputDir(OutputDir);
    }

    public boolean enabled() { return mPurpleVcf != null || PassOnly || !SpecificVariants.isEmpty(); }
    public boolean hasVariants() { return !SpecificVariants.isEmpty(); }

    public boolean hasVariantCodingImpacts() { return mPurpleVcf != null; }
    public boolean visualiserOnlyMode() { return mPurpleVcf != null; }

    public VariantContext getVariantContext(final SimpleVariant variant)
    {
        VariantAndContext variantAndContext = mReportableVariants.stream().filter(x -> x.Variant.matches(variant)).findFirst().orElse(null);
        return variantAndContext != null ? variantAndContext.Context : null;
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

        configBuilder.addConfigItem(
                PURPLE_VCF, false, "Visualiser: Purple VCF file to use for reportable variant visualsations");
    }

    public VisConfig()
    {
        OutputDir = null;
        SpecificVariants = Lists.newArrayList();
        mReportableVariants = Collections.emptyList();
        PassOnly = false;
        MaxSupportReads = 0;
        mPurpleVcf = null;
    }
}

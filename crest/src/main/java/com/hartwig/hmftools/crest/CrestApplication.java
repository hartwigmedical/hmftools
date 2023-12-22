package com.hartwig.hmftools.crest;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;

import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

// TODO arg for thresholds
// TODO output to file
public class CrestApplication
{
    private static final String RNA_SAMPLE = "rna";

    private static final Logger LOGGER = LogManager.getLogger(CrestApplication.class);

    public static void main(String[] args) throws IOException
    {
        logVersion();

        ConfigBuilder configBuilder = new ConfigBuilder();
        addConfig(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        try
        {
            String purpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
            String sampleId = configBuilder.getValue(SAMPLE);
            String rnaSample = configBuilder.getValue(RNA_SAMPLE);
            String rnaAnnotatedGermlineVcf = PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId);
            computeRnaSupportedSnpsRatio(rnaAnnotatedGermlineVcf, rnaSample);
        }
        catch(Exception e)
        {
            LOGGER.error("Crest failed: {}", e.toString());
            e.printStackTrace();
            System.exit(1);
        }
    }

    private static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addConfigItem(RNA_SAMPLE, "RNA sample ID");
        configBuilder.addConfigItem(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
    }

    public static void logVersion()
    {
        final VersionInfo version = new VersionInfo("crest.version");
        LOGGER.info("Cest version: {}", version.version());
    }

    public static double computeRnaSupportedSnpsRatio(String germlineVcf, String rnaSample) throws IOException
    {
        int supported = 0;
        var total = 0;

        try(AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(germlineVcf, new VCFCodec(), false))
        {
            for(VariantContext context : reader.iterator())
            {
                VariantContextDecorator decorator = new VariantContextDecorator(context);

                if(decorator.filter().equals("PASS") && decorator.type() == VariantType.SNP && !decorator.gene().isEmpty())
                {
                    AllelicDepth rnaDepth = decorator.allelicDepth(rnaSample);
                    if(rnaDepth.totalReadCount() >= 10)
                    {
                        total += 1;
                        if(rnaDepth.alleleReadCount() >= 1)
                        {
                            supported += 1;
                        }
                    }
                }
            }
        }
        double ratio = total > 0 ? supported * 1D / total : 0D;
        LOGGER.info("Supported: " + supported + " Total: " + total + " Fraction: " + ratio);
        return ratio;
    }
}
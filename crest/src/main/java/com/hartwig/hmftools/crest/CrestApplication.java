package com.hartwig.hmftools.crest;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.FileOutputStream;
import java.io.IOException;

import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class CrestApplication
{
    private static final Logger LOGGER = LogManager.getLogger(CrestApplication.class);

    private static final String RNA_SAMPLE = "rna_sample";
    private static final String MIN_TOTAL_READS = "min_total_reads";
    private static final String MIN_RNA_READS = "min_rna_reads";
    private static final String ACCEPTANCE_RATIO = "acceptance_ratio";
    private static final String DO_NOT_WRITE_FILE = "do_not_write_file";

    private static final String DEFAULT_MIN_TOTAL_READS = "10";
    private static final String DEFAULT_MIN_RNA_READS = "1";
    private static final String DEFAULT_ACCEPTANCE_RATIO = "0.9";

    @NotNull
    private final String purpleDir;
    @Nullable
    private final String outputDir;
    @NotNull
    private final String sampleId;
    @NotNull
    private final String sampleToCheck;

    private final int minTotalReads;
    private final int minRnaReads;
    private final double acceptanceRatio;
    private final boolean doNotWriteFile;

    public CrestApplication(@NotNull final String purpleDir, @Nullable final String outputDir,
            @NotNull final String sampleId, @NotNull final String sampleToCheck,
            final int minTotalReads, final int minRnaReads, final double acceptanceRatio,
            final boolean doNotWriteFile)
    {
        this.purpleDir = purpleDir;
        this.outputDir = outputDir;
        this.sampleId = sampleId;
        this.sampleToCheck = sampleToCheck;
        this.minTotalReads = minTotalReads;
        this.minRnaReads = minRnaReads;
        this.acceptanceRatio = acceptanceRatio;
        this.doNotWriteFile = doNotWriteFile;
    }

    public CrestApplication(@NotNull final ConfigBuilder configBuilder)
    {
        purpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        outputDir = parseOutputDir(configBuilder);
        sampleId = configBuilder.getValue(SAMPLE);
        sampleToCheck = configBuilder.getValue(RNA_SAMPLE);
        minTotalReads = configBuilder.getInteger(MIN_TOTAL_READS);
        minRnaReads = configBuilder.getInteger(MIN_RNA_READS);
        acceptanceRatio = configBuilder.getDecimal(ACCEPTANCE_RATIO);
        doNotWriteFile = configBuilder.hasFlag(DO_NOT_WRITE_FILE);
    }

    private static void registerConfig(@NotNull final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC);
        configBuilder.addConfigItem(RNA_SAMPLE, "ID of RNA sample in vcf to check, e.g. 'COLO829_RNA");
        configBuilder.addConfigItem(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addFlag(DO_NOT_WRITE_FILE, "Do not write final success or failure file");
        configBuilder.addConfigItem(MIN_TOTAL_READS, false, "Minimum total reads for a SNP to be included", DEFAULT_MIN_TOTAL_READS);
        configBuilder.addConfigItem(MIN_RNA_READS, false, "Minimum RNA reads for a SNP to be included", DEFAULT_MIN_RNA_READS);
        configBuilder.addConfigItem(ACCEPTANCE_RATIO, false, "Lower bound on fraction of RNA to total reads for check to pass", DEFAULT_ACCEPTANCE_RATIO);
        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
    }

    public static void main(String[] args) throws IOException
    {
        logVersion();

        ConfigBuilder configBuilder = new ConfigBuilder();
        registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        CrestApplication crestApplication = new CrestApplication(configBuilder);
        crestApplication.logParams();
        crestApplication.run();
    }

    public static void logVersion()
    {
        final VersionInfo version = new VersionInfo("crest.version");
        LOGGER.info("Crest version: {}", version.version());
    }

    private static boolean sampleInFile(@NotNull final String sample, @NotNull final VCFHeader header)
    {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }

    void run() throws IOException
    {
        String rnaAnnotatedGermlineVcf = PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId);
        LOGGER.info("Checking file: {}", rnaAnnotatedGermlineVcf);
        double supportRatio = computeRnaSupportRatio(rnaAnnotatedGermlineVcf);

        String outputFilename;
        if(supportRatio < acceptanceRatio)
        {
            LOGGER.error("Check failed, ratio of supported reads is below threshold");
            outputFilename = getOutputFilename(false);
        }
        else
        {
            LOGGER.info("Check succeeded");
            outputFilename = getOutputFilename(true);
        }

        if(!doNotWriteFile)
        {
            LOGGER.info("Writing file: {}", outputFilename);
            new FileOutputStream(outputFilename).close();
        }
    }

    private String getOutputFilename(boolean success)
    {
        String extension = success ? ".CrestCheckSucceeded" : ".CrestCheckFailed";
        return (outputDir == null ? "" : outputDir) + sampleId + extension;
    }

    public void logParams()
    {
        LOGGER.info("purpleDir: {}", purpleDir);
        LOGGER.info("outputDir: {}", outputDir);
        LOGGER.info("sampleId: {}", sampleId);
        LOGGER.info("sampleToCheck: {}", sampleToCheck);
        LOGGER.info("minTotalReads: {}", minTotalReads);
        LOGGER.info("minRnaReads: {}", minRnaReads);
        LOGGER.info("acceptanceRatio: {}", acceptanceRatio);
        LOGGER.info("doNotWriteFile: {}", doNotWriteFile);
    }

    public double computeRnaSupportRatio(@NotNull String vcfFilename) throws IOException
    {
        int supported = 0;
        var total = 0;

        try(AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFilename, new VCFCodec(), false))
        {
            final VCFHeader header = (VCFHeader) reader.getHeader();
            if(!sampleInFile(sampleToCheck, header))
            {
                throw new RuntimeException("Sample " + sampleToCheck + " not found in file " + vcfFilename);
            }

            for(VariantContext context : reader.iterator())
            {
                VariantContextDecorator decorator = new VariantContextDecorator(context);

                if(decorator.filter().equals("PASS") && decorator.type() == VariantType.SNP && !decorator.gene().isEmpty())
                {
                    AllelicDepth rnaDepth = decorator.allelicDepth(sampleToCheck);
                    if(rnaDepth.totalReadCount() >= minTotalReads)
                    {
                        total += 1;
                        if(rnaDepth.alleleReadCount() >= minRnaReads)
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
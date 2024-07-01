package com.hartwig.hmftools.crest;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.FileOutputStream;
import java.io.IOException;

import com.hartwig.hmftools.common.purple.PurpleCommon;
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

public class CrestAlgo
{
    private static final Logger LOGGER = LogManager.getLogger(CrestAlgo.class);

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

    public CrestAlgo(@NotNull final String purpleDir, @Nullable final String outputDir,
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

    void run() throws IOException
    {
        logVersion();
        logParams();

        String rnaAnnotatedGermlineVcf = PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId);
        LOGGER.info("Checking file: {}", rnaAnnotatedGermlineVcf);

        boolean success = crestCheck(rnaAnnotatedGermlineVcf);

        if(success)
        {
            LOGGER.info("Check succeeded");
        }
        else
        {
            LOGGER.error("Check failed, ratio of supported reads is below threshold");
        }

        if(!doNotWriteFile)
        {
            String outputFilename = getOutputFilename(success);
            LOGGER.info("Writing file: {}", outputFilename);
            new FileOutputStream(outputFilename).close();
        }
    }

    public boolean crestCheck(@NotNull String vcfFile) throws IOException
    {
        double supportRatio = computeRnaSupportRatio(vcfFile);
        return supportRatio >= acceptanceRatio;
    }

    public double computeRnaSupportRatio(@NotNull String vcfFile) throws IOException
    {
        int supported = 0;
        var total = 0;

        try(AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(vcfFile, new VCFCodec(), false))
        {
            final VCFHeader header = (VCFHeader) reader.getHeader();
            if(!sampleInFile(sampleToCheck, header))
            {
                throw new RuntimeException("Sample " + sampleToCheck + " not found in file " + vcfFile);
            }

            for(VariantContext context : reader.iterator())
            {
                VariantContextDecorator decorator = new VariantContextDecorator(context);

                if(decorator.isPass() && decorator.type() == VariantType.SNP && !decorator.gene().isEmpty())
                {
                    AllelicDepth rnaDepth = decorator.allelicDepth(sampleToCheck);
                    if(rnaDepth.TotalReadCount >= minTotalReads)
                    {
                        total += 1;
                        if(rnaDepth.AlleleReadCount >= minRnaReads)
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

    private void logParams()
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

    private static void logVersion()
    {
        final VersionInfo version = new VersionInfo("crest.version");
        LOGGER.info("Crest version: {}", version.version());
    }

    private String getOutputFilename(boolean success)
    {
        String extension = success ? ".CrestCheckSucceeded" : ".CrestCheckFailed";
        return (outputDir == null ? "" : checkAddDirSeparator(outputDir)) + sampleId + extension;
    }

    private static boolean sampleInFile(@NotNull final String sample, @NotNull final VCFHeader header)
    {
        return header.getSampleNamesInOrder().stream().anyMatch(x -> x.equals(sample));
    }
}
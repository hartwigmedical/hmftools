package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.variant.VariantContextDecorator;

public class ResultsWriter
{
    private final PurityConfig mConfig;
    private final BufferedWriter mSampleWriter;
    private final BufferedWriter mVariantWriter;

    public ResultsWriter(final PurityConfig config)
    {
        mConfig = config;
        mSampleWriter = initialiseWriter();
        mVariantWriter = config.WriteVariants ? initialiseVariantWriter() : null;
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mConfig.OutputDir + mConfig.PatientId + ".sample_summary";

            if(mConfig.OutputId != null)
                fileName += "." + mConfig.OutputId;

            fileName += ".csv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId");
            writer.write(format(",%s", CnPurityResult.header()));
            writer.write(format(",%s", SomaticVariantResult.header()));
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise variant output file: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeSampleSummary(
            final String sampleId, final CnPurityResult cnPurityResult, final SomaticVariantResult somaticVariantResult)
    {
        try
        {
            mSampleWriter.write(format("%s",sampleId));
            mSampleWriter.write(format(",%s", cnPurityResult.toCsv()));
            mSampleWriter.write(format(",%s", somaticVariantResult.toCsv()));
            mSampleWriter.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write sample output file: {}", e.toString());
        }
    }

    private BufferedWriter initialiseVariantWriter()
    {
        try
        {
            String fileName = mConfig.OutputDir + mConfig.PatientId + ".somatic_variants";

            if(mConfig.OutputId != null)
                fileName += "." + mConfig.OutputId;

            fileName += ".csv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId,Chromosome,Position,Ref,Alt,Tier,Type,RepeatCount,Mappability,SubclonalPerc,AD,DP,QualPerAD");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise variant output file: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeVariant(
            final String sampleId, final VariantContextDecorator variant,
            double subclonalLikelihood, int alleleCount, int depth, double qualPerAlleleCount)
    {
        if(mVariantWriter == null)
            return;

        try
        {
            mVariantWriter.write(format("%s,%s,%d,%s,%s",
                    sampleId, variant.chromosome(), variant.position(), variant.ref(), variant.alt()));

            mVariantWriter.write(format(",%s,%s,%d,%.3f,%.3f",
                    variant.tier(), variant.type(), variant.repeatCount(), variant.mappability(), subclonalLikelihood));

            mVariantWriter.write(format(",%d,%d,%.1f", alleleCount, depth, qualPerAlleleCount));

            mVariantWriter.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mVariantWriter);
        closeBufferedWriter(mSampleWriter);
    }

}

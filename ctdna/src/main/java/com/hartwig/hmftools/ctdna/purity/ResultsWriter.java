package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.purple.PurityContext;

public class ResultsWriter
{
    private final PurityConfig mConfig;
    private final BufferedWriter mSampleWriter;
    private final BufferedWriter mVariantWriter;
    private final BufferedWriter mCnRatioWriter;

    public ResultsWriter(final PurityConfig config)
    {
        mConfig = config;
        mSampleWriter = initialiseWriter();
        mVariantWriter = config.WriteSomatics ? initialiseVariantWriter() : null;
        mCnRatioWriter = config.WriteCnRatios ? initialiseCnRatioWriter() : null;
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mConfig.formFilename("summary");

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(mConfig.multipleSamples())
                writer.write("PatientId\t");

            writer.write("SampleId\tTumorPurity\tTumorPloidy");
            writer.write(format("\t%s", SomaticVariantResult.header()));
            writer.write(format("\t%s", CnPurityResult.header()));
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
            final String patientId, final String sampleId, final PurityContext purityContext, final CnPurityResult cnPurityResult,
            final SomaticVariantResult somaticVariantResult)
    {
        try
        {
            if(mConfig.multipleSamples())
                mSampleWriter.write(format("%s\t", patientId));

            mSampleWriter.write(format("%s\t%.2f\t%.2f", sampleId, purityContext.bestFit().purity(), purityContext.bestFit().ploidy()));
            mSampleWriter.write(format("\t%s", somaticVariantResult.toTsv()));
            mSampleWriter.write(format("\t%s", cnPurityResult.toTsv()));
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
            String fileName = mConfig.formFilename("somatic_variants");

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(mConfig.multipleSamples())
                writer.write("PatientId\t");

            writer.write("SampleId\tChromosome\tPosition\tRef\tAlt\tFilter\tTier\tType\tRepeatCount\tMappability\tSubclonalPerc\tAD\tDP\tQualPerAD");
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
            final String patientId, final String sampleId, final SomaticVariant variant, final GenotypeFragments sampleData,
            final String filter)
    {
        if(mVariantWriter == null)
            return;

        try
        {
            if(mConfig.multipleSamples())
                mVariantWriter.write(format("%s\t", patientId));

            mVariantWriter.write(format("%s\t%s\t%d\t%s\t%s\t%s",
                    sampleId, variant.Chromosome, variant.Position, variant.Ref, variant.Alt, filter));

            mVariantWriter.write(format("\t%s\t%s\t%d\t%.3f\t%.3f",
                    variant.Tier, variant.Type, variant.RepeatCount, variant.Mappability, variant.SubclonalPerc));

            mVariantWriter.write(format("\t%d\t%d\t%.1f", sampleData.AlleleCount, sampleData.Depth, sampleData.qualPerAlleleFragment()));

            mVariantWriter.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    private BufferedWriter initialiseCnRatioWriter()
    {
        try
        {
            String fileName = mConfig.formFilename("cn_segments");

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(mConfig.multipleSamples())
                writer.write("PatientId\t");

            writer.write("SampleId\tChromosome\tSegmentStart\tSegmentEnd\tCopyNumber\tGcRatioCount\tGcRatioMedian\tGcRatioMean");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise copy number segment file: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeCnSegmentData(final String patientId, final String sampleId, final CopyNumberGcData cnSegment)
    {
        if(mCnRatioWriter == null)
            return;

        // writeGcRatioData(sampleId, cnSegment); // ratios are written before they are sorted for median calcs

        try
        {
            if(mConfig.multipleSamples())
                mCnRatioWriter.write(format("%s\t", patientId));

            mCnRatioWriter.write(format("%s\t%s\t%d\t%d\t%.2f\t%d\t%.4f\t%.4f",
                    sampleId, cnSegment.Chromosome, cnSegment.SegmentStart, cnSegment.SegmentEnd, cnSegment.CopyNumber,
                    cnSegment.count(), cnSegment.median(), cnSegment.mean()));
            mCnRatioWriter.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write copy number segment file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mVariantWriter);
        closeBufferedWriter(mSampleWriter);
        closeBufferedWriter(mCnRatioWriter);
    }

    public static String formatPurityValue(double purity)
    {
        if(purity >= 0.01)
            return format("%.4f", purity);
        else
            return format("%4.3e", purity);
    }
}

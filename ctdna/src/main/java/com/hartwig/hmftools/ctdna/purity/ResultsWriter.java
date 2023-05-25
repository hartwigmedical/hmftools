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
    private final BufferedWriter mGcRatioWriter;

    public ResultsWriter(final PurityConfig config)
    {
        mConfig = config;
        mSampleWriter = initialiseWriter();
        mVariantWriter = config.WriteVariants ? initialiseVariantWriter() : null;
        mCnRatioWriter = config.WriteCnRatios ? initialiseCnRatioWriter() : null;

        mGcRatioWriter = null;
        // mGcRatioWriter = config.WriteCnRatios ? initialiseGcRatioWriter() : null;
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mConfig.formFilename("summary");

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(mConfig.multipleSamples())
                writer.write("PatientId,");

            writer.write("SampleId,TumorPurity,TumorPloidy");
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
            final String patientId, final String sampleId, final PurityContext purityContext, final CnPurityResult cnPurityResult,
            final SomaticVariantResult somaticVariantResult)
    {
        try
        {
            if(mConfig.multipleSamples())
                mSampleWriter.write(format("%s,", patientId));

            mSampleWriter.write(format("%s,%.2f,%.2f", sampleId, purityContext.bestFit().purity(), purityContext.bestFit().ploidy()));
            mSampleWriter.write(format(",%s", cnPurityResult.toTsv()));
            mSampleWriter.write(format(",%s", somaticVariantResult.toTsv()));
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
                writer.write("PatientId,");

            writer.write("SampleId,Chromosome,Position,Ref,Alt,Filter,Tier,Type,RepeatCount,Mappability,SubclonalPerc,AD,DP,QualPerAD");
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
                mVariantWriter.write(format("%s,", patientId));

            mVariantWriter.write(format("%s,%s,%d,%s,%s,%s",
                    sampleId, variant.Chromosome, variant.Position, variant.Ref, variant.Alt, filter));

            mVariantWriter.write(format(",%s,%s,%d,%.3f,%.3f",
                    variant.Tier, variant.Type, variant.RepeatCount, variant.Mappability, variant.SubclonalPerc));

            mVariantWriter.write(format(",%d,%d,%.1f", sampleData.AlleleCount, sampleData.Depth, sampleData.qualPerAlleleFragment()));

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
                writer.write("PatientId,");

            writer.write("SampleId,Chromosome,SegmentStart,SegmentEnd,CopyNumber,GcRatioCount,GcRatioMedian,GcRatioMean");
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
                mCnRatioWriter.write(format("%s,", patientId));

            mCnRatioWriter.write(format("%s,%s,%d,%d,%.2f,%d,%.4f,%.4f",
                    sampleId, cnSegment.Chromosome, cnSegment.SegmentStart, cnSegment.SegmentEnd, cnSegment.CopyNumber,
                    cnSegment.count(), cnSegment.median(), cnSegment.mean()));
            mCnRatioWriter.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write copy number segment file: {}", e.toString());
        }
    }

    private BufferedWriter initialiseGcRatioWriter()
    {
        try
        {
            String fileName = mConfig.formFilename("gc_ratio_data");

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId,Chromosome,Position,SegmentStart,SegmentEnd,CopyNumber,TumorGcRatio");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise copy number segment file: {}", e.toString());
            return null;
        }
    }

    private void writeGcRatioData(final String sampleId, final CopyNumberGcData cnSegment)
    {
        if(mGcRatioWriter == null)
            return;

        try
        {
            for(GcRatioData gcRatioData : cnSegment.ratios())
            {
                mGcRatioWriter.write(format("%s,%s,%d,%d,%d,%.4f,%.4f",
                        sampleId, cnSegment.Chromosome, gcRatioData.Position, cnSegment.SegmentStart, cnSegment.SegmentEnd,
                        cnSegment.CopyNumber, gcRatioData.TumorGcRatio));
                mGcRatioWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write copy number GC ratio file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mVariantWriter);
        closeBufferedWriter(mSampleWriter);
        closeBufferedWriter(mCnRatioWriter);
        closeBufferedWriter(mGcRatioWriter);
    }

    public static String formatPurityValue(double purity)
    {
        if(purity >= 0.01)
            return format("%.4f", purity);
        else
            return format("%4.3e", purity);
    }
}

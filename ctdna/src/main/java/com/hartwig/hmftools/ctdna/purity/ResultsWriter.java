package com.hartwig.hmftools.ctdna.purity;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.purple.PurityContext;

public class ResultsWriter
{
    private final PurityConfig mConfig;
    private final BufferedWriter mSampleWriter;
    private final BufferedWriter mVariantWriter;
    private final BufferedWriter mCnRatioWriter;

    public static final String SUMMARY_FILE_ID = "summary";
    public static final String SOMATICS_FILE_ID = "somatic_variants";
    public static final String CN_SEGMENT_FILE_ID = "cn_segments";

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
            String fileName = mConfig.formFilename(SUMMARY_FILE_ID);

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
            String fileName = mConfig.formFilename(SOMATICS_FILE_ID);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(mConfig.multipleSamples())
                sj.add("PatientId");

            sj.add("SampleId").add("Chromosome").add("Position").add("Ref").add("Alt");
            sj.add("Filter").add("Tier").add("Type").add("RepeatCount").add("Mappability").add("SubclonalPerc");
            sj.add("TumorDP").add("TumorAD");
            sj.add("SampleDP").add("SampleAD").add("SampleRefDual").add("SampleAlleleDual").add("SampleQualPerAD");
            writer.write(sj.toString());
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
            final GenotypeFragments tumorData, final String filter)
    {
        if(mVariantWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(mConfig.multipleSamples())
                sj.add(patientId);

            sj.add(sampleId).add(variant.Chromosome).add(String.valueOf(variant.Position)).add(variant.Ref).add(variant.Alt);

            sj.add(filter).add(variant.Tier.toString()).add(variant.Type.toString()).add(String.valueOf(variant.RepeatCount))
                    .add(format("%.2f", variant.Mappability)).add(format("%.2f", variant.SubclonalPerc));

            sj.add(String.valueOf(tumorData.Depth)).add(String.valueOf(tumorData.AlleleCount));
            sj.add(String.valueOf(sampleData.UmiCounts.total())).add(String.valueOf(sampleData.UmiCounts.alleleTotal()));
            sj.add(String.valueOf(sampleData.UmiCounts.RefDual)).add(String.valueOf(sampleData.UmiCounts.AlleleDual));
            sj.add(format("%.1f", sampleData.qualPerAlleleFragment()));

            mVariantWriter.write(sj.toString());

            mVariantWriter.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write output file: {}", e.toString());
            System.exit(1);
        }
    }

    private BufferedWriter initialiseCnRatioWriter()
    {
        try
        {
            String fileName = mConfig.formFilename(CN_SEGMENT_FILE_ID);

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

    public void flush()
    {
        try
        {
            mCnRatioWriter.flush();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to flush copy number segment file: {}", e.toString());
        }
    }

    public static String formatPurityValue(double purity)
    {
        if(purity >= 0.01)
            return format("%.4f", purity);
        else
            return format("%.6f", purity);
    }

    public static String formatProbabilityValue(double probability)
    {
        if(probability >= 0.03)
            return format("%.4f", probability);
        else
            return format("%4.3e", probability);
    }
}

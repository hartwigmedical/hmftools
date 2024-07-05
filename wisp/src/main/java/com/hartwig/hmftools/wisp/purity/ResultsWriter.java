package com.hartwig.hmftools.wisp.purity;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.WriteType.CN_DATA;
import static com.hartwig.hmftools.wisp.purity.WriteType.FRAG_LENGTHS;
import static com.hartwig.hmftools.wisp.purity.WriteType.LOH_DATA;
import static com.hartwig.hmftools.wisp.purity.WriteType.SOMATIC_DATA;
import static com.hartwig.hmftools.wisp.purity.loh.AmberLohCalcs.initialiseAmberLohWriter;
import static com.hartwig.hmftools.wisp.purity.cn.CopyNumberProfile.initialiseCnPlotCalcWriter;
import static com.hartwig.hmftools.wisp.purity.cn.CopyNumberProfile.initialiseCnRatioWriter;
import static com.hartwig.hmftools.wisp.purity.variant.SampleFragmentLengths.initialiseFragmentLengthWriter;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticVariants.initialiseVariantWriter;
import static com.hartwig.hmftools.wisp.purity.variant.VafPeakModel.initialiseSomaticPeakWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.wisp.purity.loh.AmberLohResult;
import com.hartwig.hmftools.wisp.purity.cn.CnPurityResult;
import com.hartwig.hmftools.wisp.purity.variant.SomaticPurityResult;

public class ResultsWriter
{
    private final PurityConfig mConfig;
    private final BufferedWriter mSampleSummaryWriter;
    private final BufferedWriter mVariantWriter;
    private final BufferedWriter mCnRatioWriter;
    private final BufferedWriter mAmberLohWriter;
    private final BufferedWriter mSomaticPeakWriter;
    private final BufferedWriter mCnPlotCalcWriter;
    private final BufferedWriter mFragLengthWriter;

    public ResultsWriter(final PurityConfig config)
    {
        mConfig = config;

        mSampleSummaryWriter = initialiseWriter();
        mVariantWriter = config.writeType(SOMATIC_DATA) ? initialiseVariantWriter(mConfig) : null;
        mCnRatioWriter = config.writeType(CN_DATA) || WriteType.plotCopyNumber(config.WriteTypes) ? initialiseCnRatioWriter(mConfig) : null;
        mAmberLohWriter = config.writeType(LOH_DATA) ? initialiseAmberLohWriter(mConfig) : null;
        mSomaticPeakWriter = WriteType.plotSomatics(config.WriteTypes) ? initialiseSomaticPeakWriter(mConfig) : null  ;
        mCnPlotCalcWriter = WriteType.plotCopyNumber(config.WriteTypes) ? initialiseCnPlotCalcWriter(mConfig) : null;
        mFragLengthWriter = config.writeType(FRAG_LENGTHS) ? initialiseFragmentLengthWriter(mConfig) : null;
    }

    public BufferedWriter getSomaticWriter() { return mVariantWriter; }
    public BufferedWriter getCnRatioWriter() { return mCnRatioWriter; }
    public BufferedWriter getSomaticPeakWriter() { return mSomaticPeakWriter; }
    public BufferedWriter getCnPlotCalcWriter() { return mCnPlotCalcWriter; }
    public BufferedWriter getAmberLohWriter() { return mAmberLohWriter; }
    public BufferedWriter getFragLengthWriter() { return mFragLengthWriter; }

    public static void addCommonHeaderFields(final StringJoiner sj, final PurityConfig config)
    {
        if(config.multiplePatients())
            sj.add("PatientId");

        if(config.multipleSamples())
            sj.add("SampleId");
    }

    public static void addCommonFields(final StringJoiner sj, final PurityConfig config, final SampleData sampleData, final String sampleId)
    {
        if(config.multiplePatients())
            sj.add(sampleData.PatientId);

        if(config.multipleSamples())
            sj.add(sampleId);
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mConfig.formFilename(FileType.SUMMARY);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            addCommonHeaderFields(sj, mConfig);

            sj.add("TumorPurity").add("TumorPloidy");

            if(mConfig.writePurityMethodData(PurityMethod.SOMATIC_VARIANT))
                sj.add(SomaticPurityResult.header());

            if(mConfig.writePurityMethodData(PurityMethod.AMBER_LOH))
                sj.add(AmberLohResult.header());

            if(mConfig.writePurityMethodData(PurityMethod.COPY_NUMBER))
                sj.add(CnPurityResult.header());

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

    public synchronized void writeSampleSummary(
            final SampleData sampleData, final String sampleId, final PurityContext purityContext, final CnPurityResult cnPurityResult,
            final SomaticPurityResult somaticPurityResult, final AmberLohResult amberLohResult)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            addCommonFields(sj, mConfig, sampleData, sampleId);

            sj.add(format("%.2f", purityContext.bestFit().purity()));
            sj.add(format("%.2f", purityContext.bestFit().ploidy()));

            if(mConfig.writePurityMethodData(PurityMethod.SOMATIC_VARIANT))
                sj.add(format("%s", somaticPurityResult.toTsv()));

            if(mConfig.writePurityMethodData(PurityMethod.AMBER_LOH))
                sj.add(format("%s", amberLohResult.toTsv()));

            if(mConfig.writePurityMethodData(PurityMethod.COPY_NUMBER))
                sj.add(format("%s", cnPurityResult.toTsv()));

            mSampleSummaryWriter.write(sj.toString());
            mSampleSummaryWriter.newLine();
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to write sample output file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mVariantWriter);
        closeBufferedWriter(mSampleSummaryWriter);
        closeBufferedWriter(mCnRatioWriter);
        closeBufferedWriter(mAmberLohWriter);
        closeBufferedWriter(mCnPlotCalcWriter);
        closeBufferedWriter(mSomaticPeakWriter);
        closeBufferedWriter(mFragLengthWriter);
    }

    public static String formatPurityValue(double purity)
    {
        if(purity >= 0.01)
            return format("%.4f", min(purity, 1.0));
        else
            return format("%.6f", min(purity, 1.0));
    }

    public static String formatProbabilityValue(double probability)
    {
        if(probability >= 0.03)
            return format("%.4f", probability);
        else
            return format("%4.3e", probability);
    }

    public static String formatDetectionResult(double estimatedPurity, double limitOfDetection)
    {
        if(limitOfDetection >= 1)
            return "NA";

        return estimatedPurity > limitOfDetection ? "TRUE" : "FALSE";
    }
}

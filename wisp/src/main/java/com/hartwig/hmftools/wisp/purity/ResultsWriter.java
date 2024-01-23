package com.hartwig.hmftools.wisp.purity;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.WriteType.CN_DATA;
import static com.hartwig.hmftools.wisp.purity.WriteType.SOMATIC_DATA;
import static com.hartwig.hmftools.wisp.purity.cn.CopyNumberProfile.initialiseCnPlotCalcWriter;
import static com.hartwig.hmftools.wisp.purity.cn.CopyNumberProfile.initialiseCnRatioWriter;
import static com.hartwig.hmftools.wisp.purity.variant.SomaticVariants.initialiseVariantWriter;
import static com.hartwig.hmftools.wisp.purity.variant.VafPeakModel.initialiseSomaticPeakWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.wisp.purity.cn.CnPurityResult;
import com.hartwig.hmftools.wisp.purity.variant.SomaticPurityResult;

public class ResultsWriter
{
    private final PurityConfig mConfig;
    private final BufferedWriter mSampleSummaryWriter;
    private final BufferedWriter mVariantWriter;
    private final BufferedWriter mCnRatioWriter;
    private final BufferedWriter mSomaticPeakWriter;
    private final BufferedWriter mCnPlotCalcWriter;

    public ResultsWriter(final PurityConfig config)
    {
        mConfig = config;

        mSampleSummaryWriter = initialiseWriter();
        mVariantWriter = config.writeType(SOMATIC_DATA) ? initialiseVariantWriter(mConfig) : null;
        mCnRatioWriter = config.writeType(CN_DATA) || WriteType.plotCopyNumber(config.WriteTypes) ? initialiseCnRatioWriter(mConfig) : null;
        mSomaticPeakWriter = WriteType.plotSomatics(config.WriteTypes) ? initialiseSomaticPeakWriter(mConfig) : null  ;
        mCnPlotCalcWriter = WriteType.plotCopyNumber(config.WriteTypes) ? initialiseCnPlotCalcWriter(mConfig) : null;
    }

    public BufferedWriter getSomaticWriter() { return mVariantWriter; }
    public BufferedWriter getCnRatioWriter() { return mCnRatioWriter; }
    public BufferedWriter getSomaticPeakWriter() { return mSomaticPeakWriter; }
    public BufferedWriter getCnPlotCalcWriter() { return mCnPlotCalcWriter; }

    public static void addCommonHeaderFields(final StringJoiner sj, final PurityConfig config)
    {
        if(config.multiplePatients())
            sj.add("PatientId");

        if(config.multipleSamples())
            sj.add("SampleId");

        if(config.hasBatchControls())
            sj.add("BatchControl");
    }

    public static void addCommonFields(final StringJoiner sj, final PurityConfig config, final SampleData sampleData, final String sampleId)
    {
        if(config.multiplePatients())
            sj.add(sampleData.PatientId);

        if(config.multipleSamples())
            sj.add(sampleId);

        if(config.hasBatchControls())
            sj.add(String.valueOf(sampleData.isBatchControl()));
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

            if(mConfig.SummaryMethodOnlyOutput && mConfig.PurityMethods.contains(PurityMethod.SOMATIC))
                sj.add(SomaticPurityResult.header());

            if(mConfig.SummaryMethodOnlyOutput && mConfig.PurityMethods.contains(PurityMethod.COPY_NUMBER))
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
            final SomaticPurityResult somaticPurityResult)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            addCommonFields(sj, mConfig, sampleData, sampleId);

            sj.add(format("%.2f", purityContext.bestFit().purity()));
            sj.add(format("%.2f", purityContext.bestFit().ploidy()));

            if(mConfig.SummaryMethodOnlyOutput && mConfig.PurityMethods.contains(PurityMethod.SOMATIC))
                sj.add(format("%s", somaticPurityResult.toTsv()));

            if(mConfig.SummaryMethodOnlyOutput && mConfig.PurityMethods.contains(PurityMethod.COPY_NUMBER))
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
        closeBufferedWriter(mCnPlotCalcWriter);
        closeBufferedWriter(mSomaticPeakWriter);
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
}

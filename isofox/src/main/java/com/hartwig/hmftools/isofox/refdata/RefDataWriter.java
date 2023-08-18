package com.hartwig.hmftools.isofox.refdata;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.VectorUtils.sumVector;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.expression.ExpectedRatesCommon.EXP_COUNT_LENGTH_HEADER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.isofox.adjusts.FragmentSize;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.expression.CategoryCountsData;

public class RefDataWriter
{
    private final RefDataConfig mConfig;
    private BufferedWriter mExpRateWriter;
    private BufferedWriter mGcRatioWriter;

    public RefDataWriter(final RefDataConfig config)
    {
        mConfig = config;

        if(mConfig.GenerateExpectedCounts)
        {
            mExpRateWriter = initExpectedCountsWriter();
        }

        if(mConfig.GenerateGcRatios)
        {
            mGcRatioWriter = initGcRatioWriter();
        }
    }

    public BufferedWriter getExpRatesWriter() { return mExpRateWriter;}
    public BufferedWriter getReadGcRatioWriter() { return mGcRatioWriter; }

    public BufferedWriter initExpectedCountsWriter()
    {
        try
        {
            String outputFileName = String.format("%sread_%d_exp_counts.%s.csv",
                    mConfig.OutputDir, mConfig.ReadLength, mConfig.RefGenVersion.identifier());

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("GeneSetId,Category");

            for(FragmentSize fragLength : mConfig.FragmentSizeData)
            {
                writer.write(String.format(",%s%d", EXP_COUNT_LENGTH_HEADER, fragLength.Length));
            }

            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected rates file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeExpectedCounts(
            final BufferedWriter writer, final String collectionId, final List<CategoryCountsData> categoryCounts)
    {
        if(writer == null)
            return;

        try
        {
            for(CategoryCountsData tcData : categoryCounts)
            {
                final double[] lengthCounts = tcData.fragmentCountsByLength();

                if(lengthCounts == null)
                    continue;

                writer.write(String.format("%s,%s,%.0f", collectionId, tcData.combinedKey(), lengthCounts[0]));

                for (int i = 1; i < lengthCounts.length; ++i)
                {
                    writer.write(String.format(",%.0f", lengthCounts[i]));
                }

                writer.newLine();
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected counts file: {}", e.toString());
        }
    }

    private BufferedWriter initGcRatioWriter()
    {
        try
        {
            String outputFileName = String.format("%sread_%d_exp_gc_ratios.%s.csv",
                    mConfig.OutputDir, mConfig.ReadLength, mConfig.RefGenVersion.identifier());

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("TransName");

            GcRatioCounts tmp = new GcRatioCounts();

            for(Double gcRatio : tmp.getRatios())
            {
                writer.write(format(",Gcr_%.2f", gcRatio));
            }

            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected GC ratio counts file: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeExpectedGcRatios(final BufferedWriter writer, final String transName, final double[] counts)
    {
        if(writer == null)
            return;

        try
        {
            writer.write(format("%s", transName));

            // convert to percentages before writing
            double frequencyTotal = sumVector(counts);

            for(Double frequency : counts)
            {
                if(frequency == 0)
                    writer.write(",0");
                else
                    writer.write(format(",%.6f", frequency/frequencyTotal));
            }

            writer.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript expected GC ratio counts file: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mExpRateWriter);
        closeBufferedWriter(mGcRatioWriter);
    }

}

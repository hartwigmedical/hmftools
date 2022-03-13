package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.compar.CommonUtils.buildComparers;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class MismatchWriter
{
    private final ComparConfig mConfig;
    private BufferedWriter mCombinedWriter;
    private final Map<Category,BufferedWriter> mCategoryWriters;

    public MismatchWriter(final ComparConfig config)
    {
        mConfig = config;
        mCombinedWriter = null;
        mCategoryWriters = Maps.newHashMap();
    }

    public boolean initialiseOutputFiles()
    {
        String filePrefix = mConfig.OutputDir;

        if(mConfig.singleSample())
            filePrefix += mConfig.SampleIds.get(0) + ".cmp.";
        else
            filePrefix += "compar_cohort.";

        if(mConfig.OutputId != null)
            filePrefix += mConfig.OutputId + ".";

        try
        {
            if(mConfig.WriteConsolidated)
            {
                String outputFile = filePrefix + "combined.csv";

                CMP_LOGGER.info("writing output results: {}", outputFile);

                mCombinedWriter = createBufferedWriter(outputFile, false);

                if(mConfig.multiSample())
                    mCombinedWriter.write("SampleId,");

                mCombinedWriter.write(Mismatch.header());
                mCombinedWriter.newLine();
            }
            else
            {
                List<ItemComparer> comparers = buildComparers(mConfig);

                for(ItemComparer comparer : comparers)
                {
                    String outputFile = filePrefix + comparer.category().toString().toLowerCase() + ".csv";

                    CMP_LOGGER.info("writing output results: {}", outputFile);

                    BufferedWriter writer = createBufferedWriter(outputFile, false);

                    if(mConfig.multiSample())
                        writer.write("SampleId,");

                    writer.write(comparer.outputHeader());
                    writer.newLine();
                    mCategoryWriters.put(comparer.category(), writer);
                }

            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to initialise Compar output files: {}", e.toString());
            return false;
        }

        return true;
    }

    public void close()
    {
        closeBufferedWriter(mCombinedWriter);
        mCategoryWriters.values().forEach(x -> closeBufferedWriter(x));
    }

    public synchronized void writeSampleMismatches(final String sampleId, final ItemComparer comparer, final List<Mismatch> mismatches)
    {
        if(mismatches.isEmpty())
            return;

        if(mCategoryWriters.isEmpty() && mCombinedWriter == null)
            return;

        try
        {
            Category category = comparer.category();
            BufferedWriter writer = mCategoryWriters.containsKey(category) ? mCategoryWriters.get(category) : mCombinedWriter;

            for(Mismatch mismatch : mismatches)
            {
                if(sampleId != null && mConfig.multiSample())
                    writer.write(String.format("%s,", sampleId));

                writer.write(comparer.mismatchOutput(mismatch));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to write sample data: {}", e.toString());
        }
    }

}

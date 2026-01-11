package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.compar.MismatchFile.loadMismatches;
import static com.hartwig.hmftools.compar.common.CommonUtils.buildComparers;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.MismatchData;

public class MismatchWriter
{
    private final ComparConfig mConfig;
    private BufferedWriter mCombinedWriter;
    private final Map<CategoryType,BufferedWriter> mCategoryWriters;

    private final Map<String,List<MismatchData>> mExpectedMismatches;

    public MismatchWriter(final ComparConfig config)
    {
        mConfig = config;
        mCombinedWriter = null;
        mCategoryWriters = Maps.newHashMap();

        String sampleId = mConfig.SampleIds.size() == 1 ? mConfig.SampleIds.get(0) : null;
        mExpectedMismatches = config.ExpectedMismatchFile != null ?
                loadMismatches(config.ExpectedMismatchFile, sampleId) : Collections.emptyMap();
    }

    public boolean initialiseOutputFiles()
    {
        String filePrefix = mConfig.OutputDir;

        if(mConfig.singleSample())
            filePrefix += mConfig.SampleIds.get(0) + ".cmp";
        else
            filePrefix += "compar_cohort";

        if(mConfig.OutputId != null)
            filePrefix += "." + mConfig.OutputId;

        try
        {
            if(mConfig.WriteDetailed)
            {
                List<ItemComparer> comparers = buildComparers(mConfig);

                for(ItemComparer comparer : comparers)
                {
                    String detailedFile = filePrefix + "." + comparer.category().toString().toLowerCase() + TSV_EXTENSION;

                    CMP_LOGGER.info("writing output results: {}", detailedFile);

                    BufferedWriter writer = createBufferedWriter(detailedFile, false);

                    writer.write(MismatchFile.commonHeader(mConfig.multiSample(), false));

                    final List<String> compareFields = comparer.comparedFieldNames();

                    for(String field : compareFields)
                    {
                        writer.write(String.format("\tRef%s\tNew%s", field, field));
                    }

                    writer.newLine();
                    mCategoryWriters.put(comparer.category(), writer);
                }
            }
            else
            {
                String outputFile = filePrefix + TSV_EXTENSION;

                CMP_LOGGER.info("writing output results: {}", outputFile);

                mCombinedWriter = createBufferedWriter(outputFile, false);

                mCombinedWriter.write(MismatchFile.header(mConfig.multiSample()));
                mCombinedWriter.newLine();
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

        checkRemoveExpectedMismatches(sampleId, mismatches);

        if(mCategoryWriters.isEmpty() && mCombinedWriter == null)
            return;

        try
        {
            CategoryType category = comparer.category();

            boolean hasSpecificWriter = mCategoryWriters.containsKey(category);
            BufferedWriter writer = hasSpecificWriter ? mCategoryWriters.get(category) : mCombinedWriter;

            for(Mismatch mismatch : mismatches)
            {
                if(sampleId != null && mConfig.multiSample())
                    writer.write(String.format("%s\t", sampleId));

                writer.write(MismatchFile.toTsv(mismatch, hasSpecificWriter, comparer.comparedFieldNames()));
                writer.newLine();
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.error("failed to write sample data: {}", e.toString());
        }
    }

    private void checkRemoveExpectedMismatches(final String sampleId, final List<Mismatch> mismatches)
    {
        List<MismatchData> expectedMismatches = mExpectedMismatches.get(sampleId);

        if(expectedMismatches == null || expectedMismatches.isEmpty())
            return;

        int index = 0;

        while(index < mismatches.size())
        {
            Mismatch mismatch = mismatches.get(index);
            boolean matched = false;

            for(int i = 0; i < expectedMismatches.size(); ++i)
            {
                MismatchData expectedMismatch = expectedMismatches.get(i);

                if(expectedMismatch.matches(mismatch))
                {
                    expectedMismatches.remove(i);
                    matched = true;
                    break;
                }
            }

            if(matched)
            {
                mismatches.remove(index);
            }
            else
            {
                ++index;
            }
        }
    }
}

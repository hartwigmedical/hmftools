package com.hartwig.hmftools.cup.prep;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_ZIP_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_CATEGORY;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_KEY;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_SOURCE;
import static com.hartwig.hmftools.cup.prep.DataItem.FLD_VALUE;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cuppa.CategoryType;

import org.jetbrains.annotations.Nullable;

public class SampleDataWriter
{
    private final PrepConfig mConfig;

    private final BufferedWriter mWriter;
    private final Map<CategoryType,BufferedWriter> mCategoryWriters;

    public SampleDataWriter(final PrepConfig config)
    {
        mConfig = config;

        if(mConfig.WriteByCategory)
        {
            mWriter = null;
            mCategoryWriters = Maps.newHashMap();

            for(CategoryType categoryType : mConfig.Categories)
            {
                BufferedWriter writer = initialiseWriter(categoryType);

                if(writer != null)
                    mCategoryWriters.put(categoryType, writer);
            }
        }
        else
        {
            mWriter = initialiseWriter(null);
            mCategoryWriters = null;
        }
    }

    public boolean isValid() { return mWriter != null || !mCategoryWriters.isEmpty(); }

    private BufferedWriter initialiseWriter(@Nullable final CategoryType categoryType)
    {
        try
        {
            String filename = mConfig.OutputDir;

            if(mConfig.isSingleSample())
            {
                filename += mConfig.SampleIds.get(0) + ".cuppa_data";
            }
            else
            {
                filename += "cuppa_cohort_data";

                if(mConfig.OutputId != null)
                {
                    filename += "." + mConfig.OutputId;
                }

                if(categoryType != null)
                {
                    filename += "." + categoryType.toString().toLowerCase();
                }
            }

            filename += TSV_ZIP_EXTENSION;

            CUP_LOGGER.info("writing {} data to {}", categoryType != null ? categoryType.toString() : "sample", filename);

            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_SOURCE);
            sj.add(FLD_CATEGORY);
            sj.add(FLD_KEY);

            if(mConfig.isSingleSample())
            {
                sj.add(FLD_VALUE);
            }
            else
            {
                for(String sampleId : mConfig.SampleIds)
                {
                    sj.add(sampleId);
                }
            }

            writer.write(sj.toString());

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to initialise writer: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeSampleData(final String sampleId, final CategoryType categoryType, final List<DataItem> dataItems)
    {
        if(mWriter == null)
            return;

        try
        {
            if(mConfig.isSingleSample())
            {
                String sourceStr = DataSource.fromCategory(categoryType).toString();

                for(DataItem dataItem : dataItems)
                {
                    mWriter.write(format("%s\t%s\t%s\t%s\t", sourceStr, dataItem.Type, dataItem.Key, dataItem.Value));
                    mWriter.newLine();
                }
            }
            else
            {
                // TODO
            }

        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write ref sample SNV counts: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mWriter);

        if(mCategoryWriters != null)
            mCategoryWriters.values().forEach(x -> closeBufferedWriter(x));
    }
}

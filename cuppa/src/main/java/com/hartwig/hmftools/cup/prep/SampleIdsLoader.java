package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.RNA_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

import org.jetbrains.annotations.Nullable;

public class SampleIdsLoader
{
    private final List<String> mSampleIds = new ArrayList<>();
    private final List<String> mRnaSampleIds = new ArrayList<>();

    public static final String COL_SAMPLE_ID = "SampleId";
    public static final String COL_RNA_SAMPLE_ID = "RnaSampleId";
    public static final String NO_RNA_SAMPLE_ID = "NO_RNA";

    public SampleIdsLoader loadFromConfig(ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(SAMPLE_ID_FILE) && !configBuilder.hasValue(SAMPLE) && !configBuilder.hasValue(RNA_SAMPLE_ID))
        {
            CUP_LOGGER.error("Sample ID(s) must be provided to 1) -{} or to 2) -{} and/or -{}", SAMPLE_ID_FILE, SAMPLE, RNA_SAMPLE_ID);
            System.exit(1);
        }

        if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            return loadFromFile(configBuilder.getValue(SAMPLE_ID_FILE));
        }
        else
        {
            return loadFromStrings(configBuilder.getValue(SAMPLE), configBuilder.getValue(RNA_SAMPLE_ID));
        }
    }

    private SampleIdsLoader loadFromStrings(String sampleId, @Nullable String rnaSampleId)
    {
        mSampleIds.add(sampleId);

        if(rnaSampleId == null)
            rnaSampleId = sampleId;

        mRnaSampleIds.add(rnaSampleId);

        return this;
    }

    @VisibleForTesting
    public SampleIdsLoader loadFromFile(String sampleIdsFile)
    {
        try
        {
            List<String> lines = FileWriterUtils.readLines(sampleIdsFile);
            return loadFromLines(lines);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("Failed to load sample IDs file: ", e);
            System.exit(1);
            return null;
        }
    }

    private SampleIdsLoader loadFromLines(List<String> lines)
    {
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        Integer sampleIdIndex = fieldsIndexMap.get(COL_SAMPLE_ID);
        Integer rnaSampleIdIndex = fieldsIndexMap.get(COL_RNA_SAMPLE_ID);
        lines.remove(0);

        if(lines.isEmpty())
        {
            CUP_LOGGER.error("No sample IDs found in sample IDs file");
            System.exit(1);
        }

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            String sampleId = values[sampleIdIndex];
            mSampleIds.add(sampleId);

            String rnaSampleId = (rnaSampleIdIndex != null)
                    ? values[rnaSampleIdIndex]
                    : "";

            if(rnaSampleId.isEmpty())
                rnaSampleId = sampleId;

            mRnaSampleIds.add(rnaSampleId);
        }

        CUP_LOGGER.debug("Loaded {} sample IDs from file", mSampleIds.size());

        return this;
    }

    public List<String> sampleIds() { return mSampleIds; }
    public List<String> rnaSampleIds() { return mRnaSampleIds; }
}

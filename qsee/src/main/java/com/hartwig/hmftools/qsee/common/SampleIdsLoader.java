package com.hartwig.hmftools.qsee.common;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.qsee.common.QseeConstants.QC_LOGGER;

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
    private final List<String> mTumorIds = new ArrayList<>();
    private final List<String> mReferenceIds = new ArrayList<>();

    private static final String COL_TUMOR_ID = "TumorId";
    private static final String COL_REFERENCE_ID = "ReferenceId";

    private static final String SAMPLE_ID_DELIM = ",";

    public SampleIdsLoader fromConfig(ConfigBuilder configBuilder)
    {
        if(!configBuilder.hasValue(SAMPLE_ID_FILE) && !configBuilder.hasValue(TUMOR) && !configBuilder.hasValue(REFERENCE))
        {
            QC_LOGGER.error("Sample IDs must be provided to 1) -{} or to 2) -{} and/or -{}", SAMPLE_ID_FILE, TUMOR, REFERENCE);
            System.exit(1);
        }

        if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            return fromFile(configBuilder.getValue(SAMPLE_ID_FILE));
        }
        else
        {
            return fromStrings(configBuilder.getValue(TUMOR), configBuilder.getValue(REFERENCE));
        }
    }

    @VisibleForTesting
    SampleIdsLoader fromStrings(String tumorIdsString, @Nullable String referenceIdsString)
    {
        mTumorIds.addAll(List.of(tumorIdsString.split(SAMPLE_ID_DELIM)));

        if(referenceIdsString != null)
            mReferenceIds.addAll(List.of(referenceIdsString.split(SAMPLE_ID_DELIM)));

        checkMatchingSampleCounts();
        return this;
    }

    private SampleIdsLoader fromFile(String sampleIdsFile)
    {
        try
        {
            List<String> lines = FileWriterUtils.readLines(sampleIdsFile);
            return fromLines(lines);
        }
        catch(IOException e)
        {
            QC_LOGGER.error("Failed to load sample IDs file:", e);
            System.exit(1);
            return null;
        }
    }

    @VisibleForTesting
    SampleIdsLoader fromLines(List<String> lines)
    {
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        Integer tumorIdIndex = fieldsIndexMap.get(COL_TUMOR_ID);
        Integer referenceIdIndex = fieldsIndexMap.get(COL_REFERENCE_ID);
        lines.remove(0);

        if(lines.isEmpty())
        {
            QC_LOGGER.error("No sample IDs found in sample IDs file");
            System.exit(1);
        }

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            String tumorId = values[tumorIdIndex];
            if(!tumorId.isEmpty())
                mTumorIds.add(tumorId);

            if(referenceIdIndex != null)
            {
                String referenceId = values[referenceIdIndex];
                if(!referenceId.isEmpty())
                    mReferenceIds.add(referenceId);
            }
        }

        checkMatchingSampleCounts();
        QC_LOGGER.debug("Loaded {} tumor and {} reference sample IDs from file", mTumorIds.size(), mReferenceIds.size());

        return this;
    }

    private void checkMatchingSampleCounts()
    {
        if(!mReferenceIds.isEmpty() && mReferenceIds.size() != mTumorIds.size())
        {
            throw new IllegalStateException("If reference samples provided, no. of reference samples must match no. of tumor samples");
        }
    }

    public List<String> tumorIds() { return mTumorIds; }
    public List<String> referenceIds() { return mReferenceIds; }
}

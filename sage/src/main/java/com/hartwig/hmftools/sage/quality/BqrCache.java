package com.hartwig.hmftools.sage.quality;

import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.redux.BqrFile;
import com.hartwig.hmftools.common.redux.BqrRecord;
import com.hartwig.hmftools.sage.SageConfig;

public class BqrCache
{
    private final SageConfig mConfig;
    private final List<String> mTumorIds;

    private final Map<String,BqrRecordMap> mSampleRecalibrationMap;
    private boolean mIsValid;

    private byte mMaxRawQual;

    public BqrCache(final SageConfig config, final List<String> tumorIds)
    {
        mConfig = config;
        mTumorIds = tumorIds;

        mSampleRecalibrationMap = Maps.newHashMap();
        mIsValid = true;
        mMaxRawQual = 0;

        if(mConfig.SkipBqr)
        {
            buildEmptyRecalibrations();
        }
        else
        {
            loadBqrFiles();
        }
    }

    public boolean isValid(){ return mIsValid; }
    public byte maxRawQual() { return mMaxRawQual; }

    public Map<String,BqrRecordMap> getSampleRecalibrationMap() { return mSampleRecalibrationMap; }

    private void loadBqrFiles()
    {
        Map<String,String> sampleFileNames = Maps.newHashMap();
        String outputDir = mConfig.outputDir();
        mConfig.ReferenceIds.forEach(x -> sampleFileNames.put(x, BqrFile.generateFilename(outputDir, x)));
        mTumorIds.forEach(x -> sampleFileNames.put(x, BqrFile.generateFilename(outputDir, x)));

        for(Map.Entry<String,String> entry : sampleFileNames.entrySet())
        {
            String sampleId = entry.getKey();
            String filename = entry.getValue();

            List<BqrRecord> bqrRecords = BqrFile.read(filename);

            if(bqrRecords == null)
            {
                mIsValid = false;
                return;
            }

            mMaxRawQual = (byte)bqrRecords.stream().mapToInt(x -> x.Key.Quality).max().orElse(0);

            SG_LOGGER.info("loaded sample({}) {} base quality recalibration records from {}", sampleId, bqrRecords.size(), filename);
            mSampleRecalibrationMap.put(sampleId, new BqrRecordMap(bqrRecords));
        }
    }

    private void buildEmptyRecalibrations()
    {
        for(String sample : mConfig.ReferenceIds)
        {
            mSampleRecalibrationMap.put(sample, new BqrRecordMap(Collections.emptyList()));
        }

        for(String sample : mTumorIds)
        {
            mSampleRecalibrationMap.put(sample, new BqrRecordMap(Collections.emptyList()));
        }
    }
}

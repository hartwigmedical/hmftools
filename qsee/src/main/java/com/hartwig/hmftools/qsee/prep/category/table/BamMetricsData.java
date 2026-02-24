package com.hartwig.hmftools.qsee.prep.category.table;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

import org.jetbrains.annotations.Nullable;

public class BamMetricsData
{
    private final CommonPrepConfig mConfig;
    private final String mSampleId;
    private final SampleType mSampleType;

    @Nullable private final BamMetricSummary mBamMetricSummary;
    @Nullable private final BamMetricCoverage mBamMetricCoverage;
    @Nullable private final BamFlagStats mBamFlagStats;

    private final List<String> mMissingInputPaths;

    public BamMetricsData(CommonPrepConfig config, String sampleId, SampleType sampleType)
    {
        mConfig = config;
        mSampleId = sampleId;
        mSampleType = sampleType;

        mMissingInputPaths = new ArrayList<>();

        mBamMetricSummary = loadBamMetricSummary();
        mBamMetricCoverage = loadBamMetricCoverage();
        mBamFlagStats = loadBamFlagStats();
    }

    @VisibleForTesting
    public BamMetricsData(@Nullable BamMetricSummary summary, BamMetricCoverage coverage, BamFlagStats flagStats, List<String> missingInputPaths)
    {
        mConfig = null;
        mSampleId = null;
        mSampleType = null;

        mMissingInputPaths = missingInputPaths;

        mBamMetricSummary = summary;
        mBamMetricCoverage = coverage;
        mBamFlagStats = flagStats;
    }

    public BamMetricSummary bamMetricSummary() { return mBamMetricSummary; }
    public BamMetricCoverage bamMetricCoverage() { return mBamMetricCoverage; }
    public BamFlagStats bamFlagStats() { return mBamFlagStats; }

    public List<String> missingInputPaths() { return mMissingInputPaths; }

    private BamMetricSummary loadBamMetricSummary()
    {
        String baseDir = mConfig.getBamMetricsDir(mSampleId, mSampleType);
        String filePath = BamMetricSummary.generateFilename(baseDir, mSampleId);

        try
        {
            return BamMetricSummary.read(filePath);
        }
        catch(IOException e)
        {
            mMissingInputPaths.add(filePath);
            return null;
        }
    }

    private BamMetricCoverage loadBamMetricCoverage()
    {
        String baseDir = mConfig.getBamMetricsDir(mSampleId, mSampleType);
        String filePath = BamMetricCoverage.generateFilename(baseDir, mSampleId);

        try
        {
            return BamMetricCoverage.read(filePath);
        }
        catch(IOException e)
        {
            mMissingInputPaths.add(filePath);
            return null;
        }
    }

    private BamFlagStats loadBamFlagStats()
    {
        String baseDir = mConfig.getBamMetricsDir(mSampleId, mSampleType);
        String filePath = BamFlagStats.generateFilename(baseDir, mSampleId);

        try
        {
            return BamFlagStats.read(filePath);
        }
        catch(IOException e)
        {
            mMissingInputPaths.add(filePath);
            return null;
        }
    }

    public String formMissingInputsString()
    {
        if(mMissingInputPaths.isEmpty())
            return "";

        StringJoiner toolsMissingInput = new StringJoiner(", ");
        for(String path : missingInputPaths())
        {
            String basename = new File(path).getName();
            toolsMissingInput.add(basename);
        }

        return toolsMissingInput.toString();
    }
}

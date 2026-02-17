package com.hartwig.hmftools.qsee.prep.category.table;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.metrics.BamFlagStats;
import com.hartwig.hmftools.common.metrics.BamMetricCoverage;
import com.hartwig.hmftools.common.metrics.BamMetricSummary;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

import org.jetbrains.annotations.Nullable;

public class SummaryTableInputData
{
    private final CommonPrepConfig mConfig;
    private final String mSampleId;
    private final SampleType mSampleType;

    @Nullable private final PurityContext mPurityContext;
    @Nullable private final BamMetricSummary mBamMetricSummary;
    @Nullable private final BamMetricCoverage mBamMetricCoverage;
    @Nullable private final BamFlagStats mBamFlagStats;

    private final List<String> mMissingInputPaths;

    public SummaryTableInputData(CommonPrepConfig config, String sampleId, SampleType sampleType)
    {
        mConfig = config;
        mSampleId = sampleId;
        mSampleType = sampleType;

        mMissingInputPaths = new ArrayList<>();

        mPurityContext = loadPurplePurity();
        mBamMetricSummary = loadBamMetricSummary();
        mBamMetricCoverage = loadBamMetricCoverage();
        mBamFlagStats = loadBamFlagStats();
    }

    public PurityContext purityContext() { return mPurityContext; }
    public BamMetricSummary bamMetricSummary() { return mBamMetricSummary; }
    public BamMetricCoverage bamMetricCoverage() { return mBamMetricCoverage; }
    public BamFlagStats bamFlagStats() { return mBamFlagStats; }

    public List<String> missingInputPaths() { return mMissingInputPaths; }

    private PurityContext loadPurplePurity()
    {
        boolean hasPurpleData = mSampleType == SampleType.TUMOR;
        if(!hasPurpleData)
            return null;

        String baseDir = mConfig.getPurpleDir(mSampleId);
        String purityFile = PurplePurity.generateFilename(baseDir, mSampleId);
        String qcFile = PurpleQCFile.generateFilename(baseDir, mSampleId);

        try
        {
            return PurityContextFile.readWithQC(qcFile, purityFile);
        }
        catch(IOException e)
        {
            mMissingInputPaths.add(purityFile);
            mMissingInputPaths.add(qcFile);
            return null;
        }
    }

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
}

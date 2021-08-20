package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.formSamplePath;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcPercentilePrevalence;
import static com.hartwig.hmftools.cup.common.CupConstants.UNDEFINED_PERC_MAX_MULTIPLE;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PERCENTILE;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadRefPercentileData;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromCohortFile;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromDatabase;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromFile;
import static com.hartwig.hmftools.cup.svs.SvDataType.LINE;
import static com.hartwig.hmftools.cup.svs.SvDataType.MAX_COMPLEX_SIZE;
import static com.hartwig.hmftools.cup.svs.SvDataType.SIMPLE_DUP_32B_200B;
import static com.hartwig.hmftools.cup.svs.SvDataType.TELOMERIC_SGL;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.linx.LinxCluster;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

public class SvClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,SvData> mSampleSvData;

    private final Map<SvDataType,Map<String,double[]>> mRefSvTypePercentiles;

    private boolean mIsValid;

    public SvClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleSvData = Maps.newHashMap();
        mRefSvTypePercentiles = Maps.newHashMap();

        mIsValid = true;

        if(mConfig.RefSvPercFile.isEmpty())
            return;

        mIsValid &= loadRefPercentileData(mConfig.RefSvPercFile, mRefSvTypePercentiles);

        CUP_LOGGER.info("loaded SV ref data from file({})", config.RefSvPercFile);

        mIsValid &= loadCohortSvData();
    }

    private static boolean isReportableType(final SvDataType type)
    {
        return (type == LINE || type == TELOMERIC_SGL || type == SIMPLE_DUP_32B_200B || type == MAX_COMPLEX_SIZE);
    }

    public CategoryType categoryType() { return SV; }
    public boolean isValid() { return mIsValid; }
    public void close() {}

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(!mIsValid || mRefSvTypePercentiles.isEmpty())
            return;

        boolean loadDbData = false;

        if(mConfig.DbAccess != null)
        {
            if(!loadSvDataFromDatabase(mConfig.DbAccess, Lists.newArrayList(sample.Id), mSampleSvData))
            {
                mIsValid = false;
                return;
            }

            loadDbData = true;
        }
        else if(mSampleSvData.isEmpty())
        {
            if(!loadSampleSvData(sample.Id))
            {
                mIsValid = false;
                return;
            }
        }

        final SvData svData = mSampleSvData.get(sample.Id);

        if(svData == null)
        {
            CUP_LOGGER.error("sample({}) sv data not loaded", sample.Id);
            mIsValid = false;
            return;
        }

        for(Map.Entry<SvDataType, Map<String, double[]>> entry : mRefSvTypePercentiles.entrySet())
        {
            final SvDataType svDataType = entry.getKey();

            if(!isReportableType(svDataType))
                continue;

            double svCount = svData.getCount(svDataType);

            final Map<String, Double> cancerTypeValues = Maps.newHashMap();

            for(Map.Entry<String, double[]> cancerPercentiles : entry.getValue().entrySet())
            {
                final String cancerType = cancerPercentiles.getKey();

                if(!isKnownCancerType(cancerType))
                    continue;

                double percentile = getPercentile(cancerPercentiles.getValue(), svCount, true, UNDEFINED_PERC_MAX_MULTIPLE);
                cancerTypeValues.put(cancerType, percentile);
            }

            SampleResult result = new SampleResult(sample.Id, SV, PERCENTILE, svDataType.toString(), svCount, cancerTypeValues);
            results.add(result);
        }

        // calculate prevalence for specific SV values
        int cancerTypeCount = mSampleDataCache.RefCancerSampleData.size();
        int cancerSampleCount = sample.isRefSample() ? mSampleDataCache.getCancerSampleCount(sample.cancerType()) : 0;

        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, LINE, true));
        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, LINE, false));
        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, TELOMERIC_SGL, false));
        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, SIMPLE_DUP_32B_200B, false));
        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, MAX_COMPLEX_SIZE, false));

        if(loadDbData)
            mSampleSvData.clear();
    }

    private SampleResult calcPrevalenceResult(
            final SampleData sample, int cancerTypeCount, int cancerSampleCount,
            final SvData svData, final SvDataType type, boolean useLowThreshold)
    {
        double svValue = svData.getCount(type);

        final Map<String,Double> cancerPrevs = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, mRefSvTypePercentiles.get(type), svValue, useLowThreshold);

        final String dataType = String.format("%s_%s", type, useLowThreshold ? "LOW" : "HIGH");
        return new SampleResult(sample.Id, SV, LIKELIHOOD, dataType, svValue, cancerPrevs);
    }

    private boolean loadCohortSvData()
    {
        if(!mConfig.SampleSvFile.isEmpty() && !mConfig.SampleSvFile.contains(".vcf"))
        {
            CUP_LOGGER.info("loading cohort SV data from file({})", mConfig.SampleSvFile);

            if(!loadSvDataFromCohortFile(mConfig.SampleSvFile, mSampleSvData))
                return false;

            CUP_LOGGER.info("loaded SV data for {} samples", mSampleSvData.size());
        }

        return true;
    }

    private boolean loadSampleSvData(final String sampleId)
    {
        if(mConfig.SampleDataDir.isEmpty() || mConfig.SampleSvFile.isEmpty() || !mConfig.SampleSvFile.contains(".vcf"))
            return false;

        final String svVcfFile = formSamplePath(mConfig.SampleSvFile, sampleId);
        final String sampleDataDir = formSamplePath(mConfig.SampleDataDir, sampleId);

        final String clusterFile = LinxCluster.generateFilename(sampleDataDir, sampleId);
        return loadSvDataFromFile(sampleId, svVcfFile, clusterFile, mSampleSvData);
    }
}

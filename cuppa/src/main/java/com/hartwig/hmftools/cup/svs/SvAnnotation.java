package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcPercentilePrevalence;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PERCENTILE;
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
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;

public class SvAnnotation
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,SvData> mSampleSvData;

    private final Map<SvDataType,Map<String,double[]>> mRefSvTypePercentiles;

    private boolean mIsValid;

    public SvAnnotation(final CuppaConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleSvData = Maps.newHashMap();
        mRefSvTypePercentiles = Maps.newHashMap();

        mIsValid = true;

        if(mConfig.RefSvPercFile.isEmpty())
            return;

        mIsValid &= loadRefPercentileData(mConfig.RefSvPercFile, mRefSvTypePercentiles);
        mIsValid &= loadSampleSvData();
    }

    private static boolean isReportableType(final SvDataType type)
    {
        return (type == LINE || type == TELOMERIC_SGL || type == SIMPLE_DUP_32B_200B || type == MAX_COMPLEX_SIZE);
    }

    public boolean isValid() { return mIsValid; }

    public List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = Lists.newArrayList();

        if(!mIsValid || mRefSvTypePercentiles.isEmpty())
            return results;

        boolean loadDbData = false;
        if(mSampleSvData.isEmpty() && mConfig.DbAccess != null)
        {
            loadDbData = true;
            if(!loadSvDataFromDatabase(mConfig.DbAccess, Lists.newArrayList(sample.Id), mSampleSvData))
                return results;
        }

        final SvData svData = mSampleSvData.get(sample.Id);

        if(svData == null)
            return results;

        for(Map.Entry<SvDataType,Map<String,double[]>> entry : mRefSvTypePercentiles.entrySet())
        {
            final SvDataType svDataType = entry.getKey();

            if(!isReportableType(svDataType))
                continue;

            double svCount = svData.getCount(svDataType);

            final Map<String,Double> cancerTypeValues = Maps.newHashMap();

            for(Map.Entry<String,double[]> cancerPercentiles : entry.getValue().entrySet())
            {
                final String cancerType = cancerPercentiles.getKey();
                double percentile = getPercentile(cancerPercentiles.getValue(), svCount, true);
                cancerTypeValues.put(cancerType, percentile);
            }

            SampleResult result = new SampleResult(sample.Id, SV, PERCENTILE, svDataType.toString(), svCount, cancerTypeValues);
            results.add(result);
        }

        // calculate prevalence for specific SV values
        int cancerTypeCount = mSampleDataCache.RefCancerSampleData.size();
        int cancerSampleCount = sample.isRefSample() ? mSampleDataCache.getCancerSampleCount(sample.CancerType) : 0;

        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, LINE, true));
        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, LINE, false));
        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, TELOMERIC_SGL, false));
        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, SIMPLE_DUP_32B_200B, false));
        results.add(calcPrevalenceResult(sample, cancerTypeCount, cancerSampleCount, svData, MAX_COMPLEX_SIZE, false));

        if(loadDbData)
            mSampleSvData.clear();

        return results;
    }

    private SampleResult calcPrevalenceResult(
            final SampleData sample, int cancerTypeCount, int cancerSampleCount,
            final SvData svData, final SvDataType type, boolean useLowThreshold)
    {
        double svValue = svData.getCount(type);

        final Map<String,Double> cancerPrevs = calcPercentilePrevalence(
                sample.CancerType, cancerSampleCount, cancerTypeCount, mRefSvTypePercentiles.get(type), svValue, useLowThreshold);

        final String dataType = String.format("%s_%s", type, useLowThreshold ? "LOW" : "HIGH");
        return new SampleResult(sample.Id, SV, LIKELIHOOD, dataType, svValue, cancerPrevs);
    }

    private boolean loadSampleSvData()
    {
        if(!mConfig.SampleSvFile.isEmpty())
        {
            if(mSampleDataCache.isSingleSample())
            {
                final String sampleId = mSampleDataCache.SampleIds.get(0);

                if(mConfig.SampleSvFile.contains(sampleId))
                {
                    final String clusterFile = LinxCluster.generateFilename(mConfig.SampleDataDir, sampleId);
                    loadSvDataFromFile(sampleId, mConfig.SampleSvFile, clusterFile, mSampleSvData);
                    return true;
                }
            }

            return loadSvDataFromCohortFile(mConfig.SampleSvFile, mSampleSvData);
        }

        return true;

        /* will load DB data for each sample
        if(mConfig.DbAccess != null)
        {
            return loadSvDataFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleSvData);
        }

        return false;
        */
    }
}

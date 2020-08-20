package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
import static com.hartwig.hmftools.cup.common.ResultType.PERCENTILE;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadRefPercentileData;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromCohortFile;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromDatabase;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.ResultType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;

import org.apache.commons.compress.utils.Lists;

public class SvAnnotation
{
    private final SampleAnalyserConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,SvData> mSampleSvData;

    private final Map<SvDataType,Map<String,double[]>> mRefSvTypePercentiles;

    private boolean mIsValid;

    public SvAnnotation(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleSvData = Maps.newHashMap();
        mRefSvTypePercentiles = Maps.newHashMap();
        mIsValid = true;

        mIsValid &= loadRefPercentileData(mConfig.RefSvPercFile, mRefSvTypePercentiles);
        mIsValid &= loadSampleSvData();
    }

    public boolean isValid() { return mIsValid; }

    public List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = Lists.newArrayList();

        if(!mConfig.runCategory(SV))
            return results;

        final SvData svData = mSampleSvData.get(sample.Id);

        if(svData == null)
            return results;

        for(Map.Entry<SvDataType,Map<String,double[]>> entry : mRefSvTypePercentiles.entrySet())
        {
            final SvDataType svDataType = entry.getKey();
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

        return results;
    }

    private boolean loadSampleSvData()
    {
        if(!mConfig.SampleSvFile.isEmpty())
        {
            if(!loadSvDataFromCohortFile(mConfig.SampleSvFile, mSampleSvData))
                return false;
        }
        else if(mConfig.DbAccess != null)
        {
            if(!loadSvDataFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleSvData))
                return false;
        }

        return true;
    }
}

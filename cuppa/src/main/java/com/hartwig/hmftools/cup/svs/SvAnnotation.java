package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadRefPercentileData;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromCohortFile;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromDatabase;
import static com.hartwig.hmftools.cup.svs.SvDataType.typeIndex;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.drivers.SampleDriverData;
import com.hartwig.hmftools.cup.sample.SampleTraitType;
import com.hartwig.hmftools.cup.sample.SampleTraitsData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.compress.utils.Lists;

public class SvAnnotation
{
    private final SampleAnalyserConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,SvData> mSampleSvData;

    private final Map<SvDataType,Map<String,double[]>> mRefSvTypePercentiles;

    public SvAnnotation(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleSvData = Maps.newHashMap();
        mRefSvTypePercentiles = Maps.newHashMap();

        loadRefPercentileData(mConfig.RefSvPercFile, mRefSvTypePercentiles);
        loadSampleSvData();
    }

    public List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = Lists.newArrayList();

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

            SampleResult result = new SampleResult(sample.Id, SV, svDataType.toString(), svCount, cancerTypeValues);
            results.add(result);
        }

        return results;
    }

    private void loadSampleSvData()
    {
        if(!mConfig.SampleSvFile.isEmpty())
        {
            loadSvDataFromCohortFile(mConfig.SampleSvFile, mSampleSvData);
        }
        else if(mConfig.DbAccess != null)
        {
            loadSvDataFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleSvData);

        }
    }

}

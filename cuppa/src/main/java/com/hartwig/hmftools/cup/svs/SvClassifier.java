package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.stats.Percentiles.getPercentile;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.common.cuppa.CategoryType.SV;
import static com.hartwig.hmftools.cup.CuppaRefFiles.purpleSvFile;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcPercentilePrevalence;
import static com.hartwig.hmftools.cup.common.CupConstants.UNDEFINED_PERC_MAX_MULTIPLE;
import static com.hartwig.hmftools.common.cuppa.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.common.cuppa.ResultType.PERCENTILE;
import static com.hartwig.hmftools.cup.common.SampleData.isKnownCancerType;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadRefPercentileData;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromCohortFile;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromDatabase;
import static com.hartwig.hmftools.cup.svs.SvDataLoader.loadSvDataFromFile;
import static com.hartwig.hmftools.common.cuppa.SvDataType.LINE;
import static com.hartwig.hmftools.common.cuppa.SvDataType.MAX_COMPLEX_SIZE;
import static com.hartwig.hmftools.common.cuppa.SvDataType.SIMPLE_DUP_32B_200B;
import static com.hartwig.hmftools.common.cuppa.SvDataType.TELOMERIC_SGL;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cuppa.SvDataType;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.common.cuppa.CategoryType;
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

    public SvClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleSvData = Maps.newHashMap();
        mRefSvTypePercentiles = Maps.newHashMap();
    }

    private static boolean isReportableType(final SvDataType type)
    {
        return (type == LINE || type == TELOMERIC_SGL || type == SIMPLE_DUP_32B_200B || type == MAX_COMPLEX_SIZE);
    }

    public CategoryType categoryType() { return SV; }
    public void close() {}

    @Override
    public boolean loadData()
    {
        if(mConfig.RefSvPercFile.isEmpty())
            return false;

        if(!loadRefPercentileData(mConfig.RefSvPercFile, mRefSvTypePercentiles))
            return false;

        CUP_LOGGER.info("loaded SV ref data from file({})", mConfig.RefSvPercFile);

        return loadSvData();
    }

    private boolean loadSvData()
    {
        if(mConfig.TestRefData)
        {
            if(!mConfig.RefSampleSvFile.isEmpty())
            {
                CUP_LOGGER.info("loading ref cohort SV data from file({})", mConfig.RefSampleSvFile);

                if(!loadSvDataFromCohortFile(mConfig.RefSampleSvFile, mSampleSvData))
                    return false;

                CUP_LOGGER.info("loaded SV data for {} samples", mSampleSvData.size());
                return true;
            }
            else
            {
                CUP_LOGGER.info("missing ref cohort SV data file");
                return false;
            }
        }

        if(mConfig.DbAccess != null)
        {
            // CUP_LOGGER.info("loading sample SV data from database"); // log if multiple only?

            return loadSvDataFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleSvData);
        }

        // CUP_LOGGER.info("loading sample SV data from pipeline files");

        for(SampleData sample : mSampleDataCache.SampleDataList)
        {
            final String purpleDataDir = mConfig.getPurpleDataDir(sample.Id);
            final String svVcfFile = purpleSvFile(purpleDataDir, sample.Id);

            final String linxDataDir = mConfig.getLinxDataDir(sample.Id);
            final String clusterFile = LinxCluster.generateFilename(linxDataDir, sample.Id, false);

            if(!loadSvDataFromFile(sample.Id, svVcfFile, clusterFile, mSampleSvData))
                return false;
        }

        return true;
    }

    @Override
    public boolean processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(mRefSvTypePercentiles.isEmpty())
            return false;

        final SvData svData = mSampleSvData.get(sample.Id);

        if(svData == null)
        {
            CUP_LOGGER.error("sample({}) sv data not loaded", sample.Id);
            return false;
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

            SampleResult result = new SampleResult(sample.Id, SV, PERCENTILE, svDataType.toString(), String.valueOf(svCount), cancerTypeValues);
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

        return true;
    }

    private SampleResult calcPrevalenceResult(
            final SampleData sample, int cancerTypeCount, int cancerSampleCount,
            final SvData svData, final SvDataType type, boolean useLowThreshold)
    {
        double svValue = svData.getCount(type);

        final Map<String,Double> cancerPrevs = calcPercentilePrevalence(
                sample, cancerSampleCount, cancerTypeCount, mRefSvTypePercentiles.get(type), svValue, useLowThreshold);

        final String dataType = String.format("%s_%s", type, useLowThreshold ? "LOW" : "HIGH");
        return new SampleResult(sample.Id, SV, LIKELIHOOD, dataType, String.valueOf(svValue), cancerPrevs);
    }

}

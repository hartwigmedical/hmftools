package com.hartwig.hmftools.cup.svs;

import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.sigs.Percentiles.getPercentile;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
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
import com.hartwig.hmftools.cup.sample.SampleTraitType;
import com.hartwig.hmftools.cup.sample.SampleTraitsData;

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

        loadRefPercentileData(mConfig.RefSvPercFile);
        loadSampleSvData(mConfig.SampleSvFile);
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

            if(svCount == 0)
                continue;

            final Map<String,Double> cancerTypeValues = Maps.newHashMap();

            for(Map.Entry<String,double[]> cancerPercentiles : entry.getValue().entrySet())
            {
                final String cancerType = cancerPercentiles.getKey();
                double percentile = getPercentile(cancerPercentiles.getValue(), svCount);
                cancerTypeValues.put(cancerType, percentile);
            }

            SampleResult result = new SampleResult(sample.Id, SV, svDataType.toString(), svCount, cancerTypeValues);
            results.add(result);
        }

        return results;
    }

    private void loadSampleSvData(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            for(final String line : fileData)
            {
                SvData svData = SvData.from(fieldsIndexMap, line);
                mSampleSvData.put(svData.SampleId, svData);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read SV data file({}): {}", filename, e.toString());
        }
    }

    private void loadRefPercentileData(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            for(final String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                final String cancerType = items[0];
                final SvDataType svDataType = SvDataType.valueOf(items[1]);

                double[] percentileData = new double[PERCENTILE_COUNT];

                int startIndex = 2;

                for(int i = startIndex; i < items.length; ++i)
                {
                    double value = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = value;
                }

                Map<String,double[]> svPercData = mRefSvTypePercentiles.get(svDataType);

                if(svPercData == null)
                {
                    svPercData = Maps.newHashMap();
                    mRefSvTypePercentiles.put(svDataType, svPercData);
                }

                svPercData.put(cancerType, percentileData);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read SV perc data file({}): {}", filename, e.toString());
        }
    }

}

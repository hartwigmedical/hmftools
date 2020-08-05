package com.hartwig.hmftools.cup.sample;

import static com.hartwig.hmftools.common.sigs.Percentiles.PERCENTILE_COUNT;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.sample.SampleTraitType.GENDER;
import static com.hartwig.hmftools.cup.sample.SampleTraitType.WGD;

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

import org.apache.commons.compress.utils.Lists;

public class SampleTraits
{
    private final SampleAnalyserConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private final Map<String,SampleTraitsData> mSampleTraitsData;

    private final Map<SampleTraitType,Map<String,double[]>> mRefTraitPercentiles;
    private final Map<SampleTraitType,Map<String,Double>> mRefTraitRates;

    public SampleTraits(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleTraitsData = Maps.newHashMap();
        mRefTraitPercentiles = Maps.newHashMap();
        mRefTraitRates = Maps.newHashMap();

        loadRefPercentileData(mConfig.RefTraitPercFile);
        loadRefRateData(mConfig.RefTraitRateFile);
        loadSampleTraitsData(mConfig.SampleTraitsFile);
    }

    public List<SampleResult> processSample(final SampleData sampleData)
    {
        List<SampleResult> results = Lists.newArrayList();

        final SampleTraitsData sampleTraits = mSampleTraitsData.get(sampleData.Id);

        if(sampleTraits == null)
            return results;

        for(Map.Entry<SampleTraitType,Map<String,Double>> entry : mRefTraitRates.entrySet())
        {
            final SampleTraitType traitType = entry.getKey();
            final Map<String,Double> cancerRates = entry.getValue();

            SampleResult result = new SampleResult(sampleData.Id, SAMPLE_TRAIT, traitType.toString(), sampleTraits.getValue(traitType), cancerRates);
            results.add(result);
        }

        return results;
    }

    private void loadSampleTraitsData(final String filename)
    {
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            for(final String line : fileData)
            {
                SampleTraitsData traitsData = SampleTraitsData.from(fieldsIndexMap, line);
                mSampleTraitsData.put(traitsData.SampleId, traitsData);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read driver prevalence data file({}): {}", filename, e.toString());
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
                final SampleTraitType traitType = SampleTraitType.valueOf(items[1]);

                double[] percentileData = new double[PERCENTILE_COUNT];

                int startIndex = 2;

                for(int i = startIndex; i < items.length; ++i)
                {
                    double value = Double.parseDouble(items[i]);
                    percentileData[i - startIndex] = value;
                }

                Map<String,double[]> traitData = mRefTraitPercentiles.get(traitType);

                if(traitData == null)
                {
                    traitData = Maps.newHashMap();
                    mRefTraitPercentiles.put(traitType, traitData);
                }

                traitData.put(cancerType, percentileData);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample traits perc data file({}): {}", filename, e.toString());
        }
    }

    private void loadRefRateData(final String filename)
    {
        // CancerType,IsFemale,WGD,SampleCount,GenderFemalePerc,WGDPerc
        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            final String header = fileData.get(0);
            fileData.remove(0);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DATA_DELIM);

            int wgdRateIndex = fieldsIndexMap.get("WGDPerc");
            int genderIndex = fieldsIndexMap.get("GenderFemalePerc");

            Map<String,Double> wgdRates = Maps.newHashMap();
            Map<String,Double> genderRates = Maps.newHashMap();

            mRefTraitRates.put(WGD, wgdRates);
            mRefTraitRates.put(GENDER, genderRates);

            for(final String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM, -1);

                final String cancerType = items[0];

                double wgdRate = Double.parseDouble(items[wgdRateIndex]);
                double genderFemaleRate = Double.parseDouble(items[genderIndex]);

                wgdRates.put(cancerType, wgdRate);
                genderRates.put(cancerType, genderFemaleRate);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample traits rate data file({}): {}", filename, e.toString());
        }
    }

}

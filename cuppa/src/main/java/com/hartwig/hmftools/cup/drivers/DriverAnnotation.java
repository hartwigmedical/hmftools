package com.hartwig.hmftools.cup.drivers;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.ClassifierType.FEATURE_PREVALENCE;
import static com.hartwig.hmftools.cup.drivers.DriverDataLoader.loadDriversFromCohortFile;
import static com.hartwig.hmftools.cup.drivers.DriverDataLoader.loadDriversFromDatabase;
import static com.hartwig.hmftools.cup.drivers.DriverDataLoader.loadRefPrevalenceData;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.ClassifierType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;

public class DriverAnnotation
{
    private final SampleAnalyserConfig mConfig;
    private final Map<String,List<SampleDriverData>> mSampleDrivers;
    private final Map<String,List<DriverPrevData>> mCancerDriverPrevalence;
    private final SampleDataCache mSampleDataCache;
    private boolean mValidData;

    private final Map<String,DriverPrevCounts> mGenePrevalenceTotals;

    public DriverAnnotation(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDrivers = Maps.newHashMap();
        mCancerDriverPrevalence = Maps.newHashMap();
        mGenePrevalenceTotals = Maps.newHashMap();
        mSampleDataCache = sampleDataCache;
        mValidData = false;

        if(config.RefDriverPrevFile.isEmpty() || config.SampleDriversFile.isEmpty())
            return;

        mValidData = true;
        loadRefPrevalenceData(config.RefDriverPrevFile, mGenePrevalenceTotals, mCancerDriverPrevalence);
        formGenePrevalenceTotals();
        loadSampleDrivers();
    }

    public final List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = Lists.newArrayList();

        if(!mValidData)
            return results;

        final List<SampleDriverData> sampleDrivers = mSampleDrivers.get(sample.Id);

        if(sampleDrivers == null || sampleDrivers.isEmpty())
            return results;

        addDriverPrevalence(sample, sampleDrivers, results);
        calcCancerTypeProbability(sample, sampleDrivers, results);

        return results;
    }

    private void addDriverPrevalence(final SampleData sample, final List<SampleDriverData> sampleDrivers, final List<SampleResult> results)
    {
        for(final SampleDriverData driverData : sampleDrivers)
        {
            final Map<String,Double> cancerTypeValues = Maps.newHashMap();

            for(Map.Entry<String,List<DriverPrevData>> entry : mCancerDriverPrevalence.entrySet())
            {
                final String cancerType = entry.getKey();

                final List<DriverPrevData> driverPrevalences = entry.getValue();

                final DriverPrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.Gene.equals(driverData.Gene)).findFirst().orElse(null);

                cancerTypeValues.put(cancerType, driverPrev != null ? driverPrev.Prevalence : 0);
            }

            SampleResult result = new SampleResult(sample.Id, FEATURE, driverData.Type.toString(), driverData.Gene, cancerTypeValues);
            results.add(result);
        }
    }

    private void calcCancerTypeProbability(
            final SampleData sample, final List<SampleDriverData> sampleDrivers, final List<SampleResult> results)
    {
        // taking the set of drivers as a group, report on the combined probability for each cancer type
        final Map<String, Double> cancerProbTotals = Maps.newHashMap();
        double allCancerProbTotal = 0;

        final Set<String> genes = Sets.newHashSet();
        sampleDrivers.forEach(x -> genes.add(x.Gene));

        final String geneNames = appendStrList(Lists.newArrayList(genes), ';');

        for(Map.Entry<String, List<DriverPrevData>> entry : mCancerDriverPrevalence.entrySet())
        {
            final String cancerType = entry.getKey();
            final List<DriverPrevData> driverPrevalences = entry.getValue();

            double probabilityTotal = 1;

            for(final String gene : genes)
            {
                double maxLikelihood = sampleDrivers.stream().filter(x -> x.Gene.equals(gene)).mapToDouble(x -> x.Likelihood).max().orElse(0);
                final DriverPrevCounts genePrevTotals = mGenePrevalenceTotals.get(gene);

                final DriverPrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                double driverPrevValue = driverPrev != null ? driverPrev.Prevalence : genePrevTotals.MinPrevalence;
                probabilityTotal *= driverPrevValue * maxLikelihood / genePrevTotals.PositiveTotal;
            }

            cancerProbTotals.put(cancerType, probabilityTotal);
            allCancerProbTotal += probabilityTotal;
        }

        final Map<String,Double> cancerTypeValues = Maps.newHashMap();

        for(Map.Entry<String,Double> entry : cancerProbTotals.entrySet())
        {
            double probability = entry.getValue() / allCancerProbTotal;
            cancerTypeValues.put(entry.getKey(), probability);
        }

        SampleResult result = new SampleResult(sample.Id, CLASSIFIER, ClassifierType.displayString(FEATURE_PREVALENCE), geneNames, cancerTypeValues);
        results.add(result);
    }


    private void loadSampleDrivers()
    {
        if(!mConfig.SampleDriversFile.isEmpty())
        {
            loadDriversFromCohortFile(mConfig.SampleDriversFile, mSampleDrivers);
        }
        else if(mConfig.DbAccess != null)
        {
            loadDriversFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleDrivers);
        }
    }

    private void formGenePrevalenceTotals()
    {
        for(Map.Entry<String,DriverPrevCounts> geneEntry : mGenePrevalenceTotals.entrySet())
        {
            final String gene = geneEntry.getKey();
            final DriverPrevCounts genePrevTotals = geneEntry.getValue();

            for(Map.Entry<String,List<DriverPrevData>> cancerEntry : mCancerDriverPrevalence.entrySet())
            {
                final DriverPrevData driverPrevData = cancerEntry.getValue().stream()
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                if(driverPrevData != null)
                {
                    genePrevTotals.PositiveTotal += driverPrevData.Prevalence;
                    genePrevTotals.NegitiveTotal += 1 - driverPrevData.Prevalence;
                }
                else
                {
                    genePrevTotals.PositiveTotal += genePrevTotals.MinPrevalence;
                    genePrevTotals.NegitiveTotal += 1 - genePrevTotals.MinPrevalence;
                }
            }
        }
    }
}

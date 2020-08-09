package com.hartwig.hmftools.cup.drivers;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.Strings.appendStr;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.cup.common.CupConstants.DRIVER_LIKELIHOOD_THRESHOLD;
import static com.hartwig.hmftools.cup.drivers.DriverDataLoader.loadFromCohortFile;
import static com.hartwig.hmftools.cup.drivers.DriverDataLoader.loadFromDatabase;
import static com.hartwig.hmftools.cup.drivers.DriverDataLoader.loadRefPrevalenceData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.sample.SampleTraitsData;

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
        loadSampleDrivers(config.SampleDriversFile);
    }

    public void processCohort()
    {
        if(!mValidData)
            return;

        mSampleDataCache.SampleDataList.forEach(x -> processSample(x));
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
                        .filter(x -> x.isTypeAll())
                        .filter(x -> x.Gene.equals(driverData.Gene)).findFirst().orElse(null);

                cancerTypeValues.put(cancerType, driverPrev != null ? driverPrev.Prevalence : 0);
            }

            SampleResult result = new SampleResult(sample.Id, DRIVER, driverData.Type, driverData.Gene, cancerTypeValues);
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

        for(Map.Entry<String, List<DriverPrevData>> entry : mCancerDriverPrevalence.entrySet())
        {
            final String cancerType = entry.getKey();
            final List<DriverPrevData> driverPrevalences = entry.getValue();

            double probabilityTotal = 1;

            for(final String driverGene : genes)
            {
                final DriverPrevCounts genePrevTotals = mGenePrevalenceTotals.get(driverGene);
                double minPrevalence = min(MIN_PREV_PERC_OF_MAX * genePrevTotals.MaxPrevalence, MIN_PREVALENCE);

                final DriverPrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.isTypeAll())
                        .filter(x -> x.Gene.equals(driverGene)).findFirst().orElse(null);

                double driverPrevValue = driverPrev != null ? driverPrev.Prevalence : minPrevalence;
                probabilityTotal *= driverPrevValue / genePrevTotals.PositiveTotal;
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

        SampleResult result = new SampleResult(sample.Id, CLASSIFIER, "DRIVER", "COMBINED", cancerTypeValues);
        results.add(result);
    }


    private void loadSampleDrivers(final String filename)
    {
        if(filename != null)
        {
            loadFromCohortFile(filename, mSampleDrivers);
        }
        else if(mConfig.DbAccess != null)
        {
            loadFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleDrivers);
        }
    }

    private static final double MIN_PREV_PERC_OF_MAX = 0.1;
    private static final double MIN_PREVALENCE = 0.01;

    private void formGenePrevalenceTotals()
    {
        for(Map.Entry<String,DriverPrevCounts> geneEntry : mGenePrevalenceTotals.entrySet())
        {
            final String gene = geneEntry.getKey();
            final DriverPrevCounts genePrevTotals = geneEntry.getValue();

            double minPrevalence = min(MIN_PREV_PERC_OF_MAX * genePrevTotals.MaxPrevalence, MIN_PREVALENCE);

            for(Map.Entry<String,List<DriverPrevData>> cancerEntry : mCancerDriverPrevalence.entrySet())
            {
                final DriverPrevData driverPrevData = cancerEntry.getValue().stream()
                        .filter(x -> x.isTypeAll())
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                if(driverPrevData != null)
                {
                    genePrevTotals.PositiveTotal += driverPrevData.Prevalence;
                    genePrevTotals.NegitiveTotal += 1 - driverPrevData.Prevalence;
                }
                else
                {
                    genePrevTotals.PositiveTotal += minPrevalence;
                    genePrevTotals.NegitiveTotal += 1 - minPrevalence;
                }
            }
        }
    }
}

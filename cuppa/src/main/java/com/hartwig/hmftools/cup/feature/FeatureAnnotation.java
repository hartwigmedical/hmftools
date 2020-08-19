package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.cup.SampleAnalyserConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.ClassifierType.FEATURE_PREVALENCE;
import static com.hartwig.hmftools.cup.common.CupConstants.DRIVER_ZERO_PREVALENCE_ALLOCATION;
import static com.hartwig.hmftools.cup.common.CupConstants.NON_DRIVER_ZERO_PREVALENCE_ALLOCATION;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadDriversFromCohortFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadDriversFromDatabase;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadRefPrevalenceData;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_CHROMOSOME;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_AMP;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_DEL;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.cup.SampleAnalyserConfig;
import com.hartwig.hmftools.cup.common.ClassifierType;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;

public class FeatureAnnotation
{
    private final SampleAnalyserConfig mConfig;
    private final Map<String,List<SampleFeatureData>> mSampleFeatures;
    private final Map<String,List<FeaturePrevData>> mCancerFeaturePrevalence;
    private final SampleDataCache mSampleDataCache;
    private boolean mValidData;

    private final Map<String, FeaturePrevCounts> mGenePrevalenceTotals;

    public FeatureAnnotation(final SampleAnalyserConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleFeatures = Maps.newHashMap();
        mCancerFeaturePrevalence = Maps.newHashMap();
        mGenePrevalenceTotals = Maps.newHashMap();
        mSampleDataCache = sampleDataCache;
        mValidData = false;

        if(config.RefFeaturePrevFile.isEmpty())
            return;

        mValidData = true;
        loadRefPrevalenceData(config.RefFeaturePrevFile, mGenePrevalenceTotals, mCancerFeaturePrevalence);
        formGenePrevalenceTotals();
        loadSampleFeatures();
    }

    public final List<SampleResult> processSample(final SampleData sample)
    {
        final List<SampleResult> results = Lists.newArrayList();

        if(!mValidData)
            return results;

        final List<SampleFeatureData> sampleFeatures = mSampleFeatures.get(sample.Id);

        if(sampleFeatures == null || sampleFeatures.isEmpty())
            return results;

        if(mConfig.runCategory(FEATURE))
        {
            addDriverPrevalence(sample, sampleFeatures, results);
        }

        if(mConfig.runCategory(CLASSIFIER))
        {
            calcCancerTypeProbability(sample, sampleFeatures, results);
        }

        return results;
    }

    private void addDriverPrevalence(final SampleData sample, final List<SampleFeatureData> sampleFeatures, final List<SampleResult> results)
    {
        for(final SampleFeatureData feature : sampleFeatures)
        {
            final Map<String,Double> cancerTypeValues = Maps.newHashMap();

            for(Map.Entry<String,List<FeaturePrevData>> entry : mCancerFeaturePrevalence.entrySet())
            {
                final String cancerType = entry.getKey();

                final List<FeaturePrevData> driverPrevalences = entry.getValue();

                final FeaturePrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.Gene.equals(feature.Gene)).findFirst().orElse(null);

                cancerTypeValues.put(cancerType, driverPrev != null ? driverPrev.Prevalence : 0);
            }

            SampleResult result = new SampleResult(sample.Id, FEATURE, feature.Type.toString(), feature.Gene, cancerTypeValues);
            results.add(result);
        }
    }

    private List<SampleFeatureData> cullMultiChromosomalEvents(final List<SampleFeatureData> features, final String cancerType)
    {
        final List<SampleFeatureData> newFeatures = Lists.newArrayList(features);
        final List<SampleFeatureData> excessFeatures = Lists.newArrayList();

        for(final SampleFeatureData feature : features)
        {
            if(feature.Type != DRIVER)
                continue;

            final String driverType = feature.getExtraInfo(DRIVER_TYPE);
            if(driverType == null || (!driverType.equals(DRIVER_TYPE_AMP) && !driverType.equals(DRIVER_TYPE_DEL)))
                continue;

            final String chromosome = feature.getExtraInfo(DRIVER_CHROMOSOME);

            // find any other matches of this driver on the same chromosome
            final List<SampleFeatureData> matchingDrivers = newFeatures.stream()
                    .filter(x -> x != feature)
                    .filter(x -> x.Type == DRIVER)
                    .filter(x -> x.getExtraInfo(DRIVER_TYPE).equals(driverType))
                    .filter(x -> !excessFeatures.contains(x))
                    .filter(x -> x.getExtraInfo(DRIVER_CHROMOSOME).equals(chromosome))
                    .collect(Collectors.toList());

            if(matchingDrivers.isEmpty())
                continue;

            matchingDrivers.add(feature);

            final List<FeaturePrevData> driverPrevalences = mCancerFeaturePrevalence.get(cancerType);

            // find the highest prevalence
            SampleFeatureData topDriver = null;
            double maxPrev = 0;
            for(SampleFeatureData matchedDriver : matchingDrivers)
            {
                final FeaturePrevData driverPrevalence = driverPrevalences.stream()
                        .filter(x -> x.Gene.equals(matchedDriver.Gene)).findFirst().orElse(null);

                if(driverPrevalence != null && driverPrevalence.Prevalence > maxPrev)
                {
                    topDriver = matchedDriver;
                    maxPrev = driverPrevalence.Prevalence;
                }
            }

            if(topDriver != null)
            {
                final SampleFeatureData topFeature = topDriver;
                matchingDrivers.stream().filter(x -> x != topFeature).forEach(x -> excessFeatures.add(x));
            }
        }

        excessFeatures.forEach(x -> newFeatures.remove(x));
        return newFeatures;
    }

    private void calcCancerTypeProbability(
            final SampleData sample, final List<SampleFeatureData> allSampleFeatures, final List<SampleResult> results)
    {
        // taking the set of drivers as a group, report on the combined probability for each cancer type
        final Map<String, Double> cancerProbTotals = Maps.newHashMap();
        double allCancerProbTotal = 0;

        final Set<String> allGenes = Sets.newHashSet();
        allSampleFeatures.forEach(x -> allGenes.add(x.Gene));
        final String geneNames = appendStrList(Lists.newArrayList(allGenes), ';');

        for(Map.Entry<String, List<FeaturePrevData>> entry : mCancerFeaturePrevalence.entrySet())
        {
            final String cancerType = entry.getKey();

            if(!sample.isCandidateCancerType(cancerType))
            {
                cancerProbTotals.put(cancerType, 0.0);
                continue;
            }

            final List<FeaturePrevData> driverPrevalences = entry.getValue();

            // only count at most one AMP or DEL per chromosome to avoid the effects of a single event impacting more than 1 gene
            final List<SampleFeatureData> sampleFeatures = cullMultiChromosomalEvents(allSampleFeatures, cancerType);

            final Set<String> genes = Sets.newHashSet();
            sampleFeatures.forEach(x -> genes.add(x.Gene));

            double probabilityTotal = 1;

            for(final String gene : genes)
            {
                double maxLikelihood = sampleFeatures.stream().filter(x -> x.Gene.equals(gene)).mapToDouble(x -> x.Likelihood).max().orElse(0);
                final FeaturePrevCounts genePrevTotals = mGenePrevalenceTotals.get(gene);

                if(genePrevTotals == null)
                {
                    CUP_LOGGER.warn("sample({}) missing gene({}) prevalence data", sample.Id, gene);
                    continue;
                }

                final FeaturePrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                double driverPrevValue = driverPrev != null ? driverPrev.Prevalence : genePrevTotals.MinPrevalence;
                probabilityTotal *= pow(driverPrevValue, maxLikelihood) / genePrevTotals.PositiveTotal;
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

    private void loadSampleFeatures()
    {
        if(!mConfig.SampleFeatureFile.isEmpty())
        {
            loadDriversFromCohortFile(mConfig.SampleFeatureFile, mSampleFeatures);
        }
        else if(mConfig.DbAccess != null)
        {
            loadDriversFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleFeatures);
        }
    }

    private void formGenePrevalenceTotals()
    {
        int cancerTypeCount = mCancerFeaturePrevalence.size();
        double noDriverPrevalence = DRIVER_ZERO_PREVALENCE_ALLOCATION / cancerTypeCount;
        double noNonDriverPrevalence = NON_DRIVER_ZERO_PREVALENCE_ALLOCATION / cancerTypeCount;

        for(Map.Entry<String,FeaturePrevCounts> geneEntry : mGenePrevalenceTotals.entrySet())
        {
            final String gene = geneEntry.getKey();
            final FeaturePrevCounts genePrevTotals = geneEntry.getValue();

            int noPrevCancerCount = 0;
            boolean isDriverType = false;

            for(Map.Entry<String,List<FeaturePrevData>> cancerEntry : mCancerFeaturePrevalence.entrySet())
            {
                final FeaturePrevData featurePrevData = cancerEntry.getValue().stream()
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                if(featurePrevData != null)
                {
                    genePrevTotals.PositiveTotal += featurePrevData.Prevalence;
                    genePrevTotals.NegitiveTotal += 1 - featurePrevData.Prevalence;
                    isDriverType = featurePrevData.Type == DRIVER;
                }
                else
                {
                    ++noPrevCancerCount;
                }
            }

            if(noPrevCancerCount > 0)
            {
                double noPrevValue = isDriverType ? noDriverPrevalence : noNonDriverPrevalence;
                genePrevTotals.PositiveTotal += noPrevValue;
                genePrevTotals.MinPrevalence = noPrevValue / noPrevCancerCount;
            }
        }
    }
}

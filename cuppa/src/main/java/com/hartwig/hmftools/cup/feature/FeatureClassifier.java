package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.max;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.ClassifierType.FEATURE_PREVALENCE;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_PAN;
import static com.hartwig.hmftools.cup.common.CupConstants.DRIVER_ZERO_PREVALENCE_ALLOCATION;
import static com.hartwig.hmftools.cup.common.CupConstants.NON_DRIVER_ZERO_PREVALENCE_ALLOCATION;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PREVALENCE;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadDriversFromCohortFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromDatabase;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadRefCancerFeatureAvg;
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
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

public class FeatureClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final Map<String,List<SampleFeatureData>> mSampleFeatures;
    private final Map<String,List<FeaturePrevData>> mCancerFeaturePrevalence;
    private final SampleDataCache mSampleDataCache;
    private boolean mIsValid;

    private final Map<String,FeaturePrevCounts> mGenePrevalenceTotals;
    private final Map<String,Double> mCancerFeatureAvg;

    public FeatureClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleFeatures = Maps.newHashMap();
        mCancerFeaturePrevalence = Maps.newHashMap();
        mGenePrevalenceTotals = Maps.newHashMap();
        mCancerFeatureAvg = Maps.newHashMap();
        mSampleDataCache = sampleDataCache;
        mIsValid = true;

        if(config.RefFeaturePrevFile.isEmpty() && config.RefDriverAvgFile.isEmpty())
            return;

        mIsValid &= loadRefPrevalenceData(config.RefFeaturePrevFile, mGenePrevalenceTotals, mCancerFeaturePrevalence);
        mIsValid &= loadRefCancerFeatureAvg(config.RefDriverAvgFile, mCancerFeatureAvg);
        formGenePrevalenceTotals();
        mIsValid &= loadSampleFeatures();
    }

    public CategoryType categoryType() { return FEATURE; }
    public boolean isValid() { return mIsValid; }

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(!mIsValid || mGenePrevalenceTotals.isEmpty())
            return;

        final List<SampleFeatureData> sampleFeatures = mSampleFeatures.get(sample.Id);

        if(sampleFeatures == null || sampleFeatures.isEmpty())
            return;

        addDriverPrevalence(sample, sampleFeatures, results);

        calcCancerTypeProbability(sample, sampleFeatures, results);
    }

    private void addDriverPrevalence(final SampleData sample, final List<SampleFeatureData> sampleFeatures, final List<SampleResult> results)
    {
        final Set<String> processedFeatures = Sets.newHashSet();

        for(final SampleFeatureData feature : sampleFeatures)
        {
            if(processedFeatures.contains(feature.Name))
                continue;

            processedFeatures.add(feature.Name);

            final Map<String,Double> cancerTypeValues = Maps.newHashMap();

            for(Map.Entry<String,List<FeaturePrevData>> entry : mCancerFeaturePrevalence.entrySet())
            {
                final String cancerType = entry.getKey();

                final List<FeaturePrevData> driverPrevalences = entry.getValue();

                final FeaturePrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.Gene.equals(feature.Name)).findFirst().orElse(null);

                cancerTypeValues.put(cancerType, driverPrev != null ? driverPrev.RawPrevalence : 0);
            }

            final String featureName = feature.Likelihood == 1 ?
                    String.format("%s (1)", feature.Name) : String.format("%s (%.2f)", feature.Name, feature.Likelihood);

            SampleResult result = new SampleResult(sample.Id, FEATURE, PREVALENCE, feature.Type.toString(), featureName, cancerTypeValues);
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
                        .filter(x -> x.Gene.equals(matchedDriver.Name)).findFirst().orElse(null);

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

    private double getFeaturesPerSampleRatio(final String cancerType)
    {
        // penalises cancer types with more features (typically drivers) per sample
        Double panCancerAvg = mCancerFeatureAvg.get(CANCER_TYPE_PAN);
        Double cancerAvg = mCancerFeatureAvg.get(cancerType);

        if(panCancerAvg != null && cancerAvg != null && panCancerAvg > 0 && cancerAvg > 0)
            return panCancerAvg / cancerAvg;

        return 1;
    }

    private void calcCancerTypeProbability(
            final SampleData sample, final List<SampleFeatureData> allSampleFeatures, final List<SampleResult> results)
    {
        // taking the set of drivers as a group, report on the combined probability for each cancer type
        final Map<String, Double> cancerProbTotals = Maps.newHashMap();

        final Set<String> allGenes = Sets.newHashSet();
        allSampleFeatures.forEach(x -> allGenes.add(x.Name));
        final String geneNames = appendStrList(Lists.newArrayList(allGenes), ';');

        for(Map.Entry<String, List<FeaturePrevData>> entry : mCancerFeaturePrevalence.entrySet())
        {
            final String cancerType = entry.getKey();

            if(!sample.isCandidateCancerType(cancerType))
            {
                cancerProbTotals.put(cancerType, 0.0);
                continue;
            }

            boolean adjustMatchingCancerPrev = sample.CancerType.equals(cancerType);

            final List<FeaturePrevData> driverPrevalences = entry.getValue();

            // only count at most one AMP or DEL per chromosome to avoid the effects of a single event impacting more than 1 gene
            final List<SampleFeatureData> sampleFeatures = allSampleFeatures; // cullMultiChromosomalEvents(allSampleFeatures, cancerType);

            final Set<String> genes = Sets.newHashSet();
            sampleFeatures.forEach(x -> genes.add(x.Name));

            double probabilityTotal = 1;

            for(final String gene : genes)
            {
                double maxLikelihood = sampleFeatures.stream().filter(x -> x.Name.equals(gene)).mapToDouble(x -> x.Likelihood).max().orElse(0);
                final FeaturePrevCounts genePrevTotals = mGenePrevalenceTotals.get(gene);

                if(genePrevTotals == null)
                {
                    CUP_LOGGER.warn("sample({}) missing gene({}) prevalence data", sample.Id, gene);
                    continue;
                }

                final FeaturePrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                double genePrevTotal = genePrevTotals.PositiveTotal;
                double driverPrevValue;

                if(driverPrev != null)
                {
                    driverPrevValue = driverPrev.Prevalence;

                    if(adjustMatchingCancerPrev)
                    {
                        int cohortSize = mSampleDataCache.getCancerSampleCount(cancerType);
                        double adjustedIncidence = driverPrevValue * cohortSize - maxLikelihood;
                        double adjustedDriverPrevValue = cohortSize > 1 ? adjustedIncidence / (cohortSize - 1) : 0;
                        genePrevTotal -= driverPrevValue - adjustedDriverPrevValue;
                        driverPrevValue = adjustedDriverPrevValue;
                    }
                }
                else
                {
                    driverPrevValue = genePrevTotals.MinPrevalence;
                }

                probabilityTotal *= pow(driverPrevValue, maxLikelihood) / genePrevTotal;
            }

            probabilityTotal *= getFeaturesPerSampleRatio(cancerType);

            cancerProbTotals.put(cancerType, probabilityTotal);
        }

        SampleResult result = new SampleResult(
                sample.Id, FEATURE, LIKELIHOOD, FEATURE_PREVALENCE.toString(), geneNames, cancerProbTotals);

        results.add(result);
    }

    private boolean loadSampleFeatures()
    {
        if(!mConfig.SampleFeatureFile.isEmpty())
        {
            return loadDriversFromCohortFile(mConfig.SampleFeatureFile, mSampleFeatures);
        }
        else if(mConfig.DbAccess != null)
        {
            return loadFeaturesFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleFeatures, false);
        }

        final String sampleId = mSampleDataCache.SampleIds.get(0);
        return loadFeaturesFromFile(sampleId, mConfig.SampleDataDir, mConfig.SampleSomaticVcf, mSampleFeatures);
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

            boolean hasNoPrevCancers = false;
            boolean isDriverType = false;

            for(Map.Entry<String,List<FeaturePrevData>> cancerEntry : mCancerFeaturePrevalence.entrySet())
            {
                final FeaturePrevData featurePrevData = cancerEntry.getValue().stream()
                        .filter(x -> x.Gene.equals(gene)).findFirst().orElse(null);

                if(featurePrevData != null)
                {
                    isDriverType = featurePrevData.Type == DRIVER;

                    double noPrevValue = isDriverType ? noDriverPrevalence : noNonDriverPrevalence;
                    featurePrevData.Prevalence += noPrevValue;
                    genePrevTotals.PositiveTotal += featurePrevData.Prevalence;
                    genePrevTotals.NegitiveTotal += 1 - featurePrevData.Prevalence;
                }
                else
                {
                    hasNoPrevCancers = true;
                }
            }

            if(hasNoPrevCancers)
            {
                double noPrevValue = isDriverType ? noDriverPrevalence : noNonDriverPrevalence;
                double noPrevAllocation = noPrevValue * cancerTypeCount;

                genePrevTotals.PositiveTotal += noPrevAllocation;
                genePrevTotals.MinPrevalence = noPrevValue;
            }
        }
    }
}

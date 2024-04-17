package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.max;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.SUBSET_DELIM;
import static com.hartwig.hmftools.common.cuppa.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_PAN;
import static com.hartwig.hmftools.cup.common.CupConstants.DRIVER_ZERO_PREVALENCE_ALLOCATION_DEFAULT;
import static com.hartwig.hmftools.cup.common.CupConstants.FEATURE_DAMPEN_FACTOR_DEFAULT;
import static com.hartwig.hmftools.cup.common.CupConstants.NON_DRIVER_ZERO_PREVALENCE_ALLOCATION_DEFAULT;
import static com.hartwig.hmftools.common.cuppa.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.common.cuppa.ResultType.PREVALENCE;
import static com.hartwig.hmftools.cup.common.CupConstants.loadKnownMutations;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromCohortFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromDatabase;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadRefCancerFeatureAvg;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadRefPrevalenceData;
import static com.hartwig.hmftools.cup.feature.FeatureType.AMP;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.FeatureType.INDEL;
import static com.hartwig.hmftools.cup.feature.FeaturesCommon.convertDriverAmps;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.cuppa.ClassifierType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

@Deprecated
public class FeatureClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final Map<String,List<SampleFeatureData>> mSampleFeatures; // keyed by sampleId
    private final Map<String,List<FeaturePrevData>> mCancerFeaturePrevalence; // keyed by cancer type
    private final SampleDataCache mSampleDataCache;

    private final Map<String,FeaturePrevCounts> mFeaturePrevalenceTotals; // keyed by feature type-name
    private final Map<String,Double> mCancerFeatureAvg; // keyed by cancer type

    private final double mNonDriverZeroPrevAllocation;
    private final double mDriverZeroPrevAllocation;

    public static final String FEATURE_DAMPEN_FACTOR = "feature_dampen_factor";
    public static final String DRIVER_ZERO_PREV = "driver_zero_prev";
    public static final String NON_DRIVER_ZERO_PREV = "non_driver_zero_prev";

    public FeatureClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final ConfigBuilder configBuilder)
    {
        mConfig = config;
        mSampleFeatures = Maps.newHashMap();
        mCancerFeaturePrevalence = Maps.newHashMap();
        mFeaturePrevalenceTotals = Maps.newHashMap();
        mCancerFeatureAvg = Maps.newHashMap();
        mSampleDataCache = sampleDataCache;

        loadKnownMutations(mConfig.RefGenVersion);

        mNonDriverZeroPrevAllocation = configBuilder.getDecimal(NON_DRIVER_ZERO_PREV);
        mDriverZeroPrevAllocation = configBuilder.getDecimal(DRIVER_ZERO_PREV);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimal(FEATURE_DAMPEN_FACTOR,"Feature dampening factor", FEATURE_DAMPEN_FACTOR_DEFAULT);

        configBuilder.addDecimal(
                DRIVER_ZERO_PREV, "Driver zero prevalence allocation", DRIVER_ZERO_PREVALENCE_ALLOCATION_DEFAULT);

        configBuilder.addDecimal(
                NON_DRIVER_ZERO_PREV,"Non-driver zero prevalence allocation", NON_DRIVER_ZERO_PREVALENCE_ALLOCATION_DEFAULT);
    }

    public CategoryType categoryType() { return FEATURE; }
    public void close() {}

    @Override
    public boolean loadData()
    {
        if(mConfig.RefFeaturePrevFile.isEmpty() && mConfig.RefDriverAvgFile.isEmpty())
            return false;

        if(!loadRefPrevalenceData(mConfig.RefFeaturePrevFile, mFeaturePrevalenceTotals, mCancerFeaturePrevalence))
            return false;

        if(!loadRefCancerFeatureAvg(mConfig.RefDriverAvgFile, mCancerFeatureAvg))
            return false;

        formFeaturePrevalenceTotals();

        CUP_LOGGER.info("loaded ref data for {} features from file({})",
                mFeaturePrevalenceTotals.size(), mConfig.RefFeaturePrevFile);

        if(!loadSampleFeatures())
            return false;

        convertDriverAmps(mSampleFeatures);

        return true;
    }

    @Override
    public boolean processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(mFeaturePrevalenceTotals.isEmpty())
            return false;

        final List<SampleFeatureData> sampleFeatures = mSampleFeatures.get(sample.Id);

        if(sampleFeatures == null || sampleFeatures.isEmpty()) // some samples legitimately have no features
            return true;

        addDriverPrevalence(sample, sampleFeatures, results);

        calcCancerTypeProbability(sample, sampleFeatures, results);
        return true;
    }

    private void addDriverPrevalence(final SampleData sample, final List<SampleFeatureData> sampleFeatures, final List<SampleResult> results)
    {
        final Set<String> processedFeatures = Sets.newHashSet();

        for(final SampleFeatureData feature : sampleFeatures)
        {
            String featureTypeName = feature.typeName();

            if(processedFeatures.contains(featureTypeName))
                continue;

            processedFeatures.add(featureTypeName);

            final Map<String,Double> cancerTypeValues = Maps.newHashMap();

            for(Map.Entry<String,List<FeaturePrevData>> entry : mCancerFeaturePrevalence.entrySet())
            {
                final String cancerType = entry.getKey();

                final List<FeaturePrevData> driverPrevalences = entry.getValue();

                final FeaturePrevData driverPrev = driverPrevalences.stream()
                        .filter(x -> x.typeName().equals(featureTypeName)).findFirst().orElse(null);

                cancerTypeValues.put(cancerType, driverPrev != null ? driverPrev.RawPrevalence : 0);
            }

            // report the max likelihood if there are multiple
            double maxLikelihood = sampleFeatures.stream()
                    .filter(x -> x.typeName().equals(featureTypeName)).mapToDouble(x -> x.Likelihood).max().orElse(0);

            final String featureName = maxLikelihood == 1 ?
                    String.format("%s (1)", feature.Name) : String.format("%s (%.2f)", feature.Name, maxLikelihood);

            SampleResult result = new SampleResult(
                    sample.Id, FEATURE, PREVALENCE, feature.Type.toString(), featureName, cancerTypeValues);
            results.add(result);
        }
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
        // skip any zero-likelihood features so they aren't given a minimum probability
        final List<SampleFeatureData> zeroLikelihoods = allSampleFeatures.stream().filter(x -> x.Likelihood == 0).collect(Collectors.toList());
        zeroLikelihoods.forEach(x -> allSampleFeatures.remove(x));

        // taking the set of drivers as a group, report on the combined probability for each cancer type
        final Map<String, Double> cancerProbTotals = Maps.newHashMap();

        final Set<String> allFeatureNames = Sets.newHashSet();

        allSampleFeatures.stream().filter(x -> x.Likelihood > 0)
                .forEach(x -> allFeatureNames.add(x.Type == AMP ? x.typeName() : x.Name));

        for(Map.Entry<String, List<FeaturePrevData>> entry : mCancerFeaturePrevalence.entrySet())
        {
            final String cancerType = entry.getKey();

            if(!checkIsValidCancerType(sample, cancerType, cancerProbTotals))
                continue;

            boolean adjustMatchingCancerPrev = sample.cancerType().equals(cancerType);

            final List<FeaturePrevData> samplePrevs = entry.getValue();

            // only count at most one driver to avoid the effects of a single event impacting more than 1 gene
            final List<SampleFeatureData> sampleFeatures = allSampleFeatures;

            final Set<String> featureTypeNames = Sets.newHashSet();
            sampleFeatures.forEach(x -> featureTypeNames.add(x.typeName()));

            double probabilityTotal = 1;

            for(final String featureTypeName : featureTypeNames)
            {
                double maxLikelihood = sampleFeatures.stream().filter(x -> x.typeName().equals(featureTypeName))
                        .mapToDouble(x -> x.Likelihood).max().orElse(0);

                if(maxLikelihood == 0)
                    continue;

                final FeaturePrevCounts featPrevTotals = mFeaturePrevalenceTotals.get(featureTypeName);

                if(featPrevTotals == null)
                {
                    CUP_LOGGER.debug("sample({}) missing gene({}) prevalence data", sample.Id, featureTypeName);
                    continue;
                }

                final FeaturePrevData featurePrev = samplePrevs.stream()
                        .filter(x -> x.typeName().equals(featureTypeName)).findFirst().orElse(null);

                double featurePrevTotal = featPrevTotals.PositiveTotal;
                double featurePrevValue;

                if(featurePrev != null)
                {
                    featurePrevValue = featurePrev.Prevalence;

                    if(adjustMatchingCancerPrev)
                    {
                        int cohortSize = mSampleDataCache.getCancerSampleCount(cancerType);
                        double origFeaturePrevValue = featurePrevValue;
                        featurePrevValue -= featPrevTotals.MinPrevalence; // remove background, then add back after recalc
                        double adjustedIncidence = max(featurePrevValue * cohortSize - maxLikelihood, 0.0);
                        featurePrevValue = cohortSize > 1 ? adjustedIncidence / (cohortSize - 1) : 0;
                        featurePrevValue += featPrevTotals.MinPrevalence;
                        featurePrevTotal -= origFeaturePrevValue - featurePrevValue;
                    }
                }
                else
                {
                    featurePrevValue = featPrevTotals.MinPrevalence;
                }

                probabilityTotal *= pow(featurePrevValue, maxLikelihood) / featurePrevTotal;
            }

            probabilityTotal *= getFeaturesPerSampleRatio(cancerType);

            cancerProbTotals.put(cancerType, probabilityTotal);
        }

        StringJoiner sj = new StringJoiner(SUBSET_DELIM);
        allFeatureNames.forEach(x -> sj.add(x));

        SampleResult result = new SampleResult(
                sample.Id, FEATURE, LIKELIHOOD, ClassifierType.FEATURE.toString(), sj.toString(), cancerProbTotals);

        results.add(result);
    }

    private boolean loadSampleFeatures()
    {
        if(mConfig.TestRefData)
        {
            if(!mConfig.RefSampleFeatureFile.isEmpty())
            {
                CUP_LOGGER.info("loading ref cohort features from file({})", mConfig.RefSampleFeatureFile);

                if(!loadFeaturesFromCohortFile(mConfig.RefSampleFeatureFile, mSampleFeatures))
                    return false;

                CUP_LOGGER.info("loaded features for {} samples", mSampleFeatures.size());
                return true;
            }
            else
            {
                CUP_LOGGER.info("missing ref cohort features file");
                return false;
            }
        }

        if(mConfig.DbAccess != null)
        {
            CUP_LOGGER.info("loading sample features from database");
            return loadFeaturesFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleFeatures);
        }

        for(SampleData sample : mSampleDataCache.SampleDataList)
        {
            final String linxDataDir = mConfig.getLinxDataDir(sample.Id);
            final String purpleDataDir = mConfig.getPurpleDataDir(sample.Id);
            final String virusDataDir = mConfig.getVirusDataDir(sample.Id);

            if(!loadFeaturesFromFile(sample.Id, linxDataDir, purpleDataDir, virusDataDir, mSampleFeatures))
                break;
        }

        return true;
    }

    private void formFeaturePrevalenceTotals()
    {
        // calculate
        int cancerTypeCount = mCancerFeaturePrevalence.size();
        double noDriverPrevalence = mDriverZeroPrevAllocation / cancerTypeCount;
        double noNonDriverPrevalence = mNonDriverZeroPrevAllocation / cancerTypeCount;

        for(Map.Entry<String,FeaturePrevCounts> prevEntry : mFeaturePrevalenceTotals.entrySet())
        {
            final String featureTypeName = prevEntry.getKey();
            final FeaturePrevCounts prevTotals = prevEntry.getValue();

            boolean isDriverType = false;

            for(Map.Entry<String,List<FeaturePrevData>> cancerEntry : mCancerFeaturePrevalence.entrySet())
            {
                final FeaturePrevData featurePrevData = cancerEntry.getValue().stream()
                        .filter(x -> x.typeName().equals(featureTypeName)).findFirst().orElse(null);

                if(featurePrevData != null)
                {
                    isDriverType = (featurePrevData.Type == DRIVER || featurePrevData.Type == AMP || featurePrevData.Type == INDEL);

                    double noPrevValue = isDriverType ? noDriverPrevalence : noNonDriverPrevalence;

                    // note adding cancer prevalence prior to added background rate, since this is done for all cancer types later
                    prevTotals.PositiveTotal += featurePrevData.Prevalence;

                    // even cancer types with non-zero prevalence are boosted by the background allocation
                    featurePrevData.Prevalence += noPrevValue;
                }
            }

            // add background rate
            double noPrevValue = isDriverType ? noDriverPrevalence : noNonDriverPrevalence;
            double noPrevAllocation = noPrevValue * cancerTypeCount;

            prevTotals.PositiveTotal += noPrevAllocation;
            prevTotals.MinPrevalence = noPrevValue;
        }
    }

    @VisibleForTesting
    public void addFeaturePrevalences(
            final Map<String,List<SampleFeatureData>> sampleFeatures,
            final Map<String,List<FeaturePrevData>> cancerFeaturePrevalence)
    {
        mSampleFeatures.putAll(sampleFeatures);
        mCancerFeaturePrevalence.putAll(cancerFeaturePrevalence);

        for(Map.Entry<String,List<FeaturePrevData>> entry : cancerFeaturePrevalence.entrySet())
        {
            String cancerType = entry.getKey();
            mCancerFeatureAvg.put(cancerType, 1.0);

            for(FeaturePrevData prevData : entry.getValue())
            {
                FeaturePrevCounts prevTotals = mFeaturePrevalenceTotals.get(prevData.typeName());

                if(prevTotals == null)
                {
                    prevTotals = new FeaturePrevCounts();
                    mFeaturePrevalenceTotals.put(prevData.typeName(), prevTotals);
                }

                prevTotals.MaxPrevalence = max(prevTotals.MaxPrevalence, prevData.Prevalence);
            }
        }

        formFeaturePrevalenceTotals();
    }
}

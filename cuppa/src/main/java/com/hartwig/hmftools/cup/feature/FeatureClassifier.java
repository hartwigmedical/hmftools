package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.max;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.SUBSET_DELIM;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CupConstants.CANCER_TYPE_PAN;
import static com.hartwig.hmftools.cup.common.CupConstants.DRIVER_ZERO_PREVALENCE_ALLOCATION;
import static com.hartwig.hmftools.cup.common.CupConstants.NON_DRIVER_ZERO_PREVALENCE_ALLOCATION;
import static com.hartwig.hmftools.cup.common.ResultType.LIKELIHOOD;
import static com.hartwig.hmftools.cup.common.ResultType.PREVALENCE;
import static com.hartwig.hmftools.cup.common.SampleResult.checkIsValidCancerType;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromCohortFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromDatabase;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadFeaturesFromFile;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadRefCancerFeatureAvg;
import static com.hartwig.hmftools.cup.feature.FeatureDataLoader.loadRefPrevalenceData;
import static com.hartwig.hmftools.cup.feature.FeatureType.AMP;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.FeatureType.INDEL;
import static com.hartwig.hmftools.cup.feature.FeaturesCommon.MIN_AMP_MULTIPLE;
import static com.hartwig.hmftools.cup.feature.FeaturesCommon.RESTRICT_DRIVER_AMP_GENES;
import static com.hartwig.hmftools.cup.feature.FeaturesCommon.COMBINE_DRIVER_AMP;
import static com.hartwig.hmftools.cup.feature.FeaturesCommon.convertAndFilterDriverAmps;
import static com.hartwig.hmftools.cup.feature.FeaturesCommon.filterDriverAmps;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.cup.CuppaConfig;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.ClassifierType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class FeatureClassifier implements CuppaClassifier
{
    private final CuppaConfig mConfig;
    private final Map<String,List<SampleFeatureData>> mSampleFeatures; // keyed by sampleId
    private final Map<String,List<FeaturePrevData>> mCancerFeaturePrevalence; // keyed by cancer type
    private final SampleDataCache mSampleDataCache;
    private boolean mIsValid;

    private final Map<String,FeaturePrevCounts> mFeaturePrevalenceTotals; // keyed by feature type-name
    private final Map<String,Double> mCancerFeatureAvg; // keyed by cancer type

    private final double mNonDriverZeroPrevAllocation;

    public static final String NON_DRIVER_ZERO_PREV = "non_driver_zero_prev";

    public FeatureClassifier(final CuppaConfig config, final SampleDataCache sampleDataCache, final CommandLine cmd)
    {
        mConfig = config;
        mSampleFeatures = Maps.newHashMap();
        mCancerFeaturePrevalence = Maps.newHashMap();
        mFeaturePrevalenceTotals = Maps.newHashMap();
        mCancerFeatureAvg = Maps.newHashMap();
        mSampleDataCache = sampleDataCache;
        mIsValid = true;

        mNonDriverZeroPrevAllocation = cmd != null && cmd.hasOption(NON_DRIVER_ZERO_PREV) ?
                Double.parseDouble(cmd.getOptionValue(NON_DRIVER_ZERO_PREV)) : NON_DRIVER_ZERO_PREVALENCE_ALLOCATION;

        if(config.RefFeaturePrevFile.isEmpty() && config.RefDriverAvgFile.isEmpty())
            return;

        mIsValid &= loadRefPrevalenceData(config.RefFeaturePrevFile, mFeaturePrevalenceTotals, mCancerFeaturePrevalence);
        mIsValid &= loadRefCancerFeatureAvg(config.RefDriverAvgFile, mCancerFeatureAvg);
        formFeaturePrevalenceTotals();

        CUP_LOGGER.info("loaded ref data for {} features from file({})",
                mFeaturePrevalenceTotals.size(), config.RefFeaturePrevFile);

        mIsValid &= loadSampleFeatures();

        boolean splitAmps = cmd == null || !cmd.hasOption(COMBINE_DRIVER_AMP);
        double minAmpCnMultiple = cmd != null ? Double.parseDouble(cmd.getOptionValue(MIN_AMP_MULTIPLE, "0")) : 0;

        if(splitAmps)
        {
            boolean restrictAmpGenes = cmd == null || cmd.hasOption(RESTRICT_DRIVER_AMP_GENES);
            convertAndFilterDriverAmps(mSampleFeatures, restrictAmpGenes);

            if(minAmpCnMultiple > 0)
            {
                filterDriverAmps(mSampleFeatures, mSampleDataCache.RefSampleTraitsData, minAmpCnMultiple);
            }
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(NON_DRIVER_ZERO_PREV, true, "Non-driver zero prevalence allocation");
        FeaturesCommon.addCmdLineArgs(options);
    }

    public CategoryType categoryType() { return FEATURE; }
    public boolean isValid() { return mIsValid; }
    public void close() {}

    public void processSample(final SampleData sample, final List<SampleResult> results, final List<SampleSimilarity> similarities)
    {
        if(!mIsValid || mFeaturePrevalenceTotals.isEmpty())
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
        if(!mConfig.SampleFeatureFile.isEmpty())
        {
            CUP_LOGGER.info("loading sample features from file({})", mConfig.SampleFeatureFile);

            if(!loadFeaturesFromCohortFile(mConfig.SampleFeatureFile, mSampleFeatures))
                return false;

            CUP_LOGGER.info("loaded features for {} samples", mSampleFeatures.size());
            return true;
        }
        else if(mConfig.DbAccess != null)
        {
            CUP_LOGGER.info("loading sample features from database");
            return loadFeaturesFromDatabase(mConfig.DbAccess, mSampleDataCache.SampleIds, mSampleFeatures);
        }

        for(SampleData sample : mSampleDataCache.SampleDataList)
        {
            final String linxDataDir = mConfig.getLinxDataDir(sample.Id);
            final String purpleDataDir = mConfig.getPurpleDataDir(sample.Id);

            if(!loadFeaturesFromFile(sample.Id, linxDataDir, purpleDataDir, mSampleFeatures))
                break;
        }

        return true;
    }

    private void formFeaturePrevalenceTotals()
    {
        // calculate
        int cancerTypeCount = mCancerFeaturePrevalence.size();
        double noDriverPrevalence = DRIVER_ZERO_PREVALENCE_ALLOCATION / cancerTypeCount;
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

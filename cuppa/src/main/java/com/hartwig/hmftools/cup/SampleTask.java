package com.hartwig.hmftools.cup;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.cuppa.DataTypes.DATA_TYPE_COMBINED;
import static com.hartwig.hmftools.common.cuppa.DataTypes.DATA_TYPE_DNA_COMBINED;
import static com.hartwig.hmftools.common.cuppa.DataTypes.DATA_TYPE_RNA_COMBINED;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcCombinedClassifierScoreResult;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcCombinedFeatureResult;
import static com.hartwig.hmftools.cup.common.CupCalcs.fillMissingCancerTypeValues;
import static com.hartwig.hmftools.cup.common.CupConstants.COMBINED_DAMPEN_FACTOR;
import static com.hartwig.hmftools.cup.common.CupConstants.DNA_DAMPEN_FACTOR;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_DAMPEN_FACTOR;
import static com.hartwig.hmftools.common.cuppa.ResultType.CLASSIFIER;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

public class SampleTask implements Callable
{
    private final int mTaskId;
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;
    private final List<CuppaClassifier> mClassifiers;
    private final ResultsWriter mResultsWriter;
    private final List<SampleData> mSamples;

    public SampleTask(
            final int taskId, final CuppaConfig config, final SampleDataCache sampleDataCache,
            final List<CuppaClassifier> classifiers, final ResultsWriter resultsWriter)
    {
        mTaskId = taskId;
        mConfig = config;
        mSampleDataCache = sampleDataCache;
        mClassifiers = classifiers;
        mResultsWriter = resultsWriter;

        mSamples = Lists.newArrayList();
    }

    @Override
    public Long call() throws Exception
    {
        run();
        return (long)0;
    }

    public List<SampleData> getSamples() { return mSamples; }

    public void run() throws Exception
    {
        int sampleCount = 0;
        for(SampleData sample : mSamples)
        {
            CUP_LOGGER.debug("sample({}) running CUP analysis", sample.Id);

            if(!processSample(sample))
            {
                // CUP_LOGGER.info("{}: exiting on error", mTaskId);
                // return;
                throw new Exception(format("task(%d) exiting on error", mTaskId));
            }

            ++sampleCount;

            if((sampleCount % 100) == 0)
            {
                CUP_LOGGER.info("{}: processed {} samples", mTaskId, sampleCount);
            }
        }

        CUP_LOGGER.info("{}: task complete", mTaskId);
    }

    public boolean processSample(final SampleData sample)
    {
        final List<SampleResult> allResults = Lists.newArrayList();
        final List<SampleSimilarity> similarities = Lists.newArrayList();

        for(CuppaClassifier classifier : mClassifiers)
        {
            if(!classifier.processSample(sample, allResults, similarities))
                return false;
        }

        // combine all features into a single classifier
        SampleResult combinedFeatureResult = calcCombinedFeatureResult(sample, allResults, mSampleDataCache.SampleIds.size() > 1);

        if(combinedFeatureResult != null)
            allResults.add(combinedFeatureResult);

        boolean hasDnaCategories = mConfig.Categories.stream().anyMatch(x -> CategoryType.isDna(x));
        boolean hasRnaCategories = mConfig.Categories.stream().anyMatch(x -> CategoryType.isRna(x));

        // ensure each cancer type has a probability for the classifiers to ensure an even application of the min-probability
        final Set<String> refCancerTypes = mSampleDataCache.RefCancerSampleData.keySet().stream()
                .filter(x -> !mSampleDataCache.RefCancerMappings.containsKey(x))
                .collect(Collectors.toSet());

        allResults.stream()
                .filter(x -> x.Result == CLASSIFIER)
                .forEach(x -> fillMissingCancerTypeValues(x.CancerTypeValues, refCancerTypes));

        if(hasDnaCategories)
        {
            final List<SampleResult> dnaResults = allResults.stream()
                    .filter(x -> x.Result == CLASSIFIER && CategoryType.isDna(x.Category))
                    .collect(Collectors.toList());

            if(dnaResults.isEmpty())
            {
                hasDnaCategories = false;
            }
            else
            {
                SampleResult dnaScoreResult = calcCombinedClassifierScoreResult(sample, dnaResults, DATA_TYPE_DNA_COMBINED, DNA_DAMPEN_FACTOR);

                if(dnaScoreResult != null)
                    allResults.add(dnaScoreResult);
            }
        }

        if(hasRnaCategories)
        {
            final List<SampleResult> rnaResults = allResults.stream()
                    .filter(x -> x.Result == CLASSIFIER && CategoryType.isRna(x.Category))
                    .collect(Collectors.toList());

            if(rnaResults.isEmpty())
            {
                hasRnaCategories = false;
            }
            else
            {
                SampleResult rnaScoreResult = calcCombinedClassifierScoreResult(sample, rnaResults, DATA_TYPE_RNA_COMBINED, RNA_DAMPEN_FACTOR);

                if(rnaScoreResult != null)
                    allResults.add(rnaScoreResult);
            }
        }

        if(hasDnaCategories && hasRnaCategories)
        {
            SampleResult classifierScoreResult =
                    calcCombinedClassifierScoreResult(sample, allResults, DATA_TYPE_COMBINED, COMBINED_DAMPEN_FACTOR);

            if(classifierScoreResult != null)
                allResults.add(classifierScoreResult);
        }

        // collapse any sub-types into parent types
        if(!mSampleDataCache.RefCancerMappings.isEmpty() && !mConfig.NoSubtypeCollapse)
        {
            collapseCancerSubtypes(allResults);
        }

        mResultsWriter.writeSampleData(sample, allResults);
        mResultsWriter.writeSampleSimilarities(sample, similarities);

        return true;
    }

    private void collapseCancerSubtypes(final List<SampleResult> results)
    {
        for(SampleResult result : results)
        {
            if(result.Result != CLASSIFIER)
                continue;

            final Map<String,Double> cancerTypeResults = result.CancerTypeValues;

            for(Map.Entry<String,String> mapping : mSampleDataCache.RefCancerMappings.entrySet())
            {
                Double subtypeProb = cancerTypeResults.get(mapping.getKey());
                Double mainTypeProb = cancerTypeResults.get(mapping.getValue());

                if(subtypeProb != null && mainTypeProb != null)
                {
                    cancerTypeResults.put(mapping.getValue(), subtypeProb + mainTypeProb);
                    cancerTypeResults.remove(mapping.getKey());
                }
            }
        }
    }

}

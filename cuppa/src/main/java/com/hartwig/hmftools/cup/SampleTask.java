package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.cup.CupAnalyser.allClassifiersValid;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.ClassifierType.isDna;
import static com.hartwig.hmftools.cup.common.ClassifierType.isRna;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcCombinedClassifierScoreResult;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcCombinedFeatureResult;
import static com.hartwig.hmftools.cup.common.CupCalcs.fillMissingCancerTypeValues;
import static com.hartwig.hmftools.cup.common.CupConstants.COMBINED_DAMPEN_FACTOR;
import static com.hartwig.hmftools.cup.common.CupConstants.DNA_DAMPEN_FACTOR;
import static com.hartwig.hmftools.cup.common.CupConstants.RNA_DAMPEN_FACTOR;

import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.ClassifierType;
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
    public Long call()
    {
        run();
        return (long)0;
    }

    public List<SampleData> getSamples() { return mSamples; }

    public void run()
    {
        int sampleCount = 0;
        for(SampleData sample : mSamples)
        {
            CUP_LOGGER.debug("sample({}) running CUP analysis", sample.Id);

            processSample(sample);

            if(!allClassifiersValid(mClassifiers))
                break;

            ++sampleCount;

            if((sampleCount % 100) == 0)
            {
                CUP_LOGGER.info("{}: processed {} samples", mTaskId, sampleCount);
            }
        }

        CUP_LOGGER.info("{}: task complete", mTaskId);
    }

    public void processSample(final SampleData sample)
    {
        final List<SampleResult> allResults = Lists.newArrayList();
        final List<SampleSimilarity> similarities = Lists.newArrayList();

        for(CuppaClassifier classifier : mClassifiers)
        {
            classifier.processSample(sample, allResults, similarities);
        }

        // combine all features into a single classifier
        SampleResult combinedFeatureResult = calcCombinedFeatureResult(sample, allResults, mSampleDataCache.SampleIds.size() > 1);

        if(combinedFeatureResult != null)
            allResults.add(combinedFeatureResult);

        boolean hasDnaCategories = mConfig.Categories.stream().anyMatch(x -> CategoryType.isDna(x));
        boolean hasRnaCategories = mConfig.Categories.stream().anyMatch(x -> CategoryType.isRna(x));

        // ensure each cancer type has a probability for the classifiers to ensure an even application of the min-probability
        final Set<String> refCancerTypes = mSampleDataCache.RefCancerSampleData.keySet();

        allResults.stream()
                .filter(x -> x.Category == CLASSIFIER)
                .forEach(x -> fillMissingCancerTypeValues(x.CancerTypeValues, refCancerTypes));

        if(hasDnaCategories)
        {
            final List<SampleResult> dnaResults = allResults.stream()
                    .filter(x -> x.Category == CLASSIFIER && isDna(ClassifierType.valueOf(x.DataType)))
                    .collect(Collectors.toList());

            if(dnaResults.isEmpty())
            {
                hasDnaCategories = false;
            }
            else
            {
                SampleResult dnaScoreResult = calcCombinedClassifierScoreResult(sample, dnaResults, "DNA_COMBINED", DNA_DAMPEN_FACTOR);

                if(dnaScoreResult != null)
                    allResults.add(dnaScoreResult);
            }
        }

        if(hasRnaCategories)
        {
            final List<SampleResult> rnaResults = allResults.stream()
                    .filter(x -> x.Category == CLASSIFIER && isRna(ClassifierType.valueOf(x.DataType)))
                    .collect(Collectors.toList());

            if(rnaResults.isEmpty())
            {
                hasRnaCategories = false;
            }
            else
            {
                SampleResult rnaScoreResult = calcCombinedClassifierScoreResult(sample, rnaResults, "RNA_COMBINED", RNA_DAMPEN_FACTOR);

                if(rnaScoreResult != null)
                    allResults.add(rnaScoreResult);
            }
        }

        if(hasDnaCategories && hasRnaCategories)
        {
            SampleResult classifierScoreResult =
                    calcCombinedClassifierScoreResult(sample, allResults, "COMBINED", COMBINED_DAMPEN_FACTOR);

            if(classifierScoreResult != null)
                allResults.add(classifierScoreResult);
        }

        mResultsWriter.writeSampleData(sample, allResults);
        mResultsWriter.writeSampleSimilarities(sample, similarities);
    }

}

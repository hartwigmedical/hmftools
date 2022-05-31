package com.hartwig.hmftools.cup;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.COMBINED;
import static com.hartwig.hmftools.cup.common.CategoryType.isDna;
import static com.hartwig.hmftools.cup.common.CategoryType.isRna;
import static com.hartwig.hmftools.cup.common.ClassifierType.ALT_SJ_COHORT;
import static com.hartwig.hmftools.cup.common.ClassifierType.EXPRESSION_PAIRWISE;
import static com.hartwig.hmftools.cup.common.ClassifierType.FEATURE;
import static com.hartwig.hmftools.cup.common.ClassifierType.GENOMIC_POSITION_COHORT;
import static com.hartwig.hmftools.cup.common.ClassifierType.SNV_96_PAIRWISE;
import static com.hartwig.hmftools.cup.common.CupConstants.DATA_TYPE_COMBINED;
import static com.hartwig.hmftools.cup.common.CupConstants.DATA_TYPE_DNA_COMBINED;
import static com.hartwig.hmftools.cup.common.CupConstants.DATA_TYPE_RNA_COMBINED;
import static com.hartwig.hmftools.cup.common.ResultType.CLASSIFIER;
import static com.hartwig.hmftools.cup.utils.CompareUtils.getRankedCancerTypes;
import static com.hartwig.hmftools.cup.utils.CompareUtils.topRefResult;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

public class ResultsWriter
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;

    private BufferedWriter mDetailedWriter;
    private BufferedWriter mCondensedWriter;
    private BufferedWriter mSimilarityWriter;

    public ResultsWriter(final CuppaConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mDetailedWriter = null;
        mCondensedWriter = null;
        mSimilarityWriter = null;

        initialiseOutputFiles();
    }

    private void initialiseOutputFiles()
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            final String detailedFilename = mSampleDataCache.isSingleSample() ?
                    mConfig.OutputDir + mSampleDataCache.SpecificSample.Id + ".cup.data.csv"
                    : mConfig.formOutputFilename("SAMPLE_DATA");

            mDetailedWriter = createBufferedWriter(detailedFilename, false);

            mDetailedWriter.write(SampleResult.detailedHeader());

            mDetailedWriter.newLine();

            if(mConfig.WriteCondensed)
            {
                final String condensedFilename = mSampleDataCache.isSingleSample() ?
                        mConfig.OutputDir + mSampleDataCache.SpecificSample.Id + ".cup.results.csv"
                        : mConfig.formOutputFilename("RESULTS");

                mCondensedWriter = createBufferedWriter(condensedFilename, false);

                mCondensedWriter.write("SampleId,Platform,CancerType,CancerTypeRank,Combined1,Combined2,Combined3");
                mCondensedWriter.write(",CombinedScore,CombinedScore1,CombinedScore2,CombinedScore3");

                mCondensedWriter.write(",DnaCombinedScore,TopDnaCombined,TopDnaCombinedScore");
                mCondensedWriter.write(",RnaCombinedScore,TopRnaCombined,TopRnaCombinedScore");
                mCondensedWriter.write(",Snv96Score,TopSnv96,TopSnv96Score,GenPosScore,TopGenPos,TopGenPosScore");
                mCondensedWriter.write(",FeatureScore,TopFeature,TopFeatureScore,ExpressionScore,TopExpression,TopExpressionScore");
                mCondensedWriter.write(",AltSjScore,TopAltSj,TopAltSjScore");

                mCondensedWriter.newLine();
            }

            if(mConfig.WriteSimilarities)
            {
                final String sampleSimilarityFilename = mSampleDataCache.isSingleSample() ?
                        mConfig.OutputDir + mSampleDataCache.SpecificSample.Id + ".cup.similarities.csv"
                        : mConfig.formOutputFilename("SAMPLE_SIMILARITIES");

                mSimilarityWriter = createBufferedWriter(sampleSimilarityFilename, false);

                mSimilarityWriter.write("SampleId,CancerType,PrimaryType,PrimarySubtype,Location,SubLocation");
                mSimilarityWriter.write(",MatchType,Score,MatchSampleId,MatchCancerType");
                mSimilarityWriter.write(",MatchPrimaryType,MatchPrimarySubtype,MatchLocation,MatchSubLocation");
                mSimilarityWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write CUPPA output: {}", e.toString());
        }
    }

    public synchronized void writeSampleData(final SampleData sample, final List<SampleResult> results)
    {
        if(results.isEmpty() || mDetailedWriter == null)
            return;

        try
        {
            for(SampleResult result : results)
            {
                if(!mConfig.WriteDetailedScores && result.Result != CLASSIFIER)
                    continue;

                result.writeDetailed(mDetailedWriter);
            }

            writeCondensedSampleResults(sample, results);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample data: {}", e.toString());
        }
    }

    private void writeCondensedSampleResults(final SampleData sample, final List<SampleResult> results) throws IOException
    {
        if(mCondensedWriter == null)
            return;

        List<SampleResult> combinedResults = results.stream().filter(x -> x.Category == COMBINED).collect(Collectors.toList());

        if(combinedResults.isEmpty())
            return;

        List<SampleResult> classifierResults = results.stream().filter(x -> x.Result == CLASSIFIER).collect(Collectors.toList());

        boolean hasDna = classifierResults.stream().anyMatch(x -> isDna(x.Category));
        boolean hasRna = classifierResults.stream().anyMatch(x -> isRna(x.Category));

        String platform;
        SampleResult combinedResult;

        if(hasDna && hasRna)
        {
            combinedResult = combinedResults.stream().filter(x -> x.DataType.equals(DATA_TYPE_COMBINED)).findFirst().orElse(null);
            platform = "WGS+WTS";
        }
        else if(hasDna)
        {
            combinedResult = combinedResults.stream().filter(x -> x.DataType.equals(DATA_TYPE_DNA_COMBINED)).findFirst().orElse(null);
            platform = "WGS";
        }
        else
        {
            combinedResult = combinedResults.stream().filter(x -> x.DataType.equals(DATA_TYPE_RNA_COMBINED)).findFirst().orElse(null);
            platform = "WTS";
        }

        List<String> rankedCancerTypes = getRankedCancerTypes(combinedResult);
        int refRank = 0;
        for(; refRank < rankedCancerTypes.size(); ++refRank)
        {
            if(rankedCancerTypes.get(refRank).equals(sample.cancerMainType()))
                break;
        }

        // SampleId,RefCancerType,RefRank,Platform

        final String sampleCancerType = sample.cancerMainType();
        mCondensedWriter.write(format("%s,%s,%s,%d", sample.Id, platform, sampleCancerType, refRank + 1));

        mCondensedWriter.write(format(",%s,%s,%s,%.4f,%.4f,%.4f,%.4f",
                rankedCancerTypes.get(0), rankedCancerTypes.get(1), rankedCancerTypes.get(2),
                combinedResult.CancerTypeValues.get(sampleCancerType),
                combinedResult.CancerTypeValues.get(rankedCancerTypes.get(0)),
                combinedResult.CancerTypeValues.get(rankedCancerTypes.get(1)),
                combinedResult.CancerTypeValues.get(rankedCancerTypes.get(2))));

        List<String> classifiers = Lists.newArrayList(
                DATA_TYPE_DNA_COMBINED, DATA_TYPE_RNA_COMBINED, SNV_96_PAIRWISE.toString(), GENOMIC_POSITION_COHORT.toString(),
                FEATURE.toString(), EXPRESSION_PAIRWISE.toString(), ALT_SJ_COHORT.toString());

        for(String classifier : classifiers)
        {
            SampleResult result = classifierResults.stream().filter(x -> x.DataType.equals(classifier))
                    .findFirst().orElse(null);

            if(result != null)
            {
                String topCancerType = topRefResult(result);

                mCondensedWriter.write(format(",%.4f,%s,%.4f",
                        result.CancerTypeValues.get(sampleCancerType), topCancerType, result.CancerTypeValues.get(topCancerType)));
            }
            else
            {
                mCondensedWriter.write(",0,N/A,0");
            }
        }

        mCondensedWriter.newLine();
    }

    public synchronized void writeSampleSimilarities(final SampleData sampleData, final List<SampleSimilarity> similarities)
    {
        if(similarities.isEmpty() || mSimilarityWriter == null)
            return;

        try
        {
            for(SampleSimilarity similarity : similarities)
            {
                SampleData matchedSample = mSampleDataCache.findRefSampleData(similarity.MatchedSampleId);

                if(matchedSample == null)
                {
                    matchedSample = mSampleDataCache.findSampleData(similarity.MatchedSampleId);
                }

                mSimilarityWriter.write(String.format("%s,%s,%s,%s,%s,%s,%s,%.3f,%s",
                        sampleData.Id, sampleData.cancerType(), sampleData.PrimaryType, sampleData.PrimarySubtype,
                        sampleData.PrimaryLocation, sampleData.PrimarySubLocation,
                        similarity.MatchType, similarity.Score, similarity.MatchedSampleId));

                if(matchedSample != null)
                {
                    mSimilarityWriter.write(String.format(",%s,%s,%s,%s,%s",
                            matchedSample.cancerType(), matchedSample.PrimaryType, matchedSample.PrimarySubtype,
                            matchedSample.PrimaryLocation, matchedSample.PrimarySubLocation));
                }
                else
                {
                    mSimilarityWriter.write(",Unclassifed,,,,");
                }

                mSimilarityWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample similarity: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mDetailedWriter);
        closeBufferedWriter(mSimilarityWriter);
        closeBufferedWriter(mCondensedWriter);
    }

}

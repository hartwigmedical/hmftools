package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.common.CategoryType.isSummary;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;

public class ResultsWriter
{
    private final CuppaConfig mConfig;
    private final SampleDataCache mSampleDataCache;
    private BufferedWriter mSampleDataWriter;
    private BufferedWriter mSampleSimilarityWriter;

    public ResultsWriter(final CuppaConfig config, final SampleDataCache sampleDataCache)
    {
        mConfig = config;
        mSampleDataCache = sampleDataCache;

        mSampleDataWriter = null;
        mSampleSimilarityWriter = null;
        initialiseOutputFiles();
    }

    private void initialiseOutputFiles()
    {
        if(mConfig.OutputDir == null)
            return;

        try
        {
            final String sampleDataFilename = mSampleDataCache.isSingleSample() ?
                    mConfig.OutputDir + mSampleDataCache.SpecificSample.Id + ".cup.data.csv"
                    : mConfig.formOutputFilename("SAMPLE_DATA");

            mSampleDataWriter = createBufferedWriter(sampleDataFilename, false);

            mSampleDataWriter.write(SampleResult.csvHeader());

            mSampleDataWriter.newLine();

            if(mConfig.WriteSimilarities)
            {
                final String sampleSimilarityFilename = mSampleDataCache.isSingleSample() ?
                        mConfig.OutputDir + mSampleDataCache.SpecificSample.Id + ".cup.similarities.csv"
                        : mConfig.formOutputFilename("SAMPLE_SIMILARITIES");

                mSampleSimilarityWriter = createBufferedWriter(sampleSimilarityFilename, false);

                mSampleSimilarityWriter.write("SampleId,CancerType,PrimaryType,PrimarySubtype,Location,SubLocation");
                mSampleSimilarityWriter.write(",MatchType,Score,MatchSampleId,MatchCancerType");
                mSampleSimilarityWriter.write(",MatchPrimaryType,MatchPrimarySubtype,MatchLocation,MatchSubLocation");
                mSampleSimilarityWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write CUPPA output: {}", e.toString());
        }
    }

    public synchronized void writeSampleData(final SampleData sampleData, final List<SampleResult> results)
    {
        if(results.isEmpty() || mSampleDataWriter == null)
            return;

        try
        {
            for(SampleResult result : results)
            {
                if(!mConfig.WriteDetailedScores && !isSummary(result.Category))
                    continue;

                result.write(mSampleDataWriter);
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample data: {}", e.toString());
        }
    }

    public synchronized void writeSampleSimilarities(final SampleData sampleData, final List<SampleSimilarity> similarities)
    {
        if(similarities.isEmpty() || mSampleSimilarityWriter == null)
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

                mSampleSimilarityWriter.write(String.format("%s,%s,%s,%s,%s,%s,%s,%.3f,%s",
                        sampleData.Id, sampleData.cancerType(), sampleData.PrimaryType, sampleData.PrimarySubtype,
                        sampleData.PrimaryLocation, sampleData.PrimarySubLocation,
                        similarity.MatchType, similarity.Score, similarity.MatchedSampleId));

                if(matchedSample != null)
                {
                    mSampleSimilarityWriter.write(String.format(",%s,%s,%s,%s,%s",
                            matchedSample.cancerType(), matchedSample.PrimaryType, matchedSample.PrimarySubtype,
                            matchedSample.PrimaryLocation, matchedSample.PrimarySubLocation));
                }
                else
                {
                    mSampleSimilarityWriter.write(",Unclassifed,,,,");
                }

                mSampleSimilarityWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample similarity: {}", e.toString());
        }
    }

    public void close()
    {
        closeBufferedWriter(mSampleDataWriter);
        closeBufferedWriter(mSampleSimilarityWriter);
    }

}

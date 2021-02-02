package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.cup.CuppaConfig.LOG_DEBUG;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.SPECIFIC_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.common.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.cup.common.CategoryType.CLASSIFIER;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;
import static com.hartwig.hmftools.cup.common.ClassifierType.COMBINED;
import static com.hartwig.hmftools.cup.common.ClassifierType.isDna;
import static com.hartwig.hmftools.cup.common.ClassifierType.isRna;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcClassifierScoreResult;
import static com.hartwig.hmftools.cup.common.CupCalcs.calcCombinedFeatureResult;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.ClassifierType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.common.SampleResult;
import com.hartwig.hmftools.cup.common.SampleSimilarity;
import com.hartwig.hmftools.cup.feature.FeatureClassifier;
import com.hartwig.hmftools.cup.rna.AltSjClassifier;
import com.hartwig.hmftools.cup.rna.GeneExpressionClassifier;
import com.hartwig.hmftools.cup.sample.SampleTraitClassifier;
import com.hartwig.hmftools.cup.somatics.SomaticClassifier;
import com.hartwig.hmftools.cup.svs.SvClassifier;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class CupAnalyser
{
    private final CuppaConfig mConfig;

    private final SampleDataCache mSampleDataCache;

    private final List<CuppaClassifier> mClassifiers;

    private BufferedWriter mSampleDataWriter;
    private BufferedWriter mSampleSimilarityWriter;

    public CupAnalyser(final CommandLine cmd)
    {
        mConfig = new CuppaConfig(cmd);

        mSampleDataCache = new SampleDataCache();

        loadSampleData(cmd);

        mClassifiers = Lists.newArrayList();

        if(mConfig.runClassifier(SNV))
            mClassifiers.add(new SomaticClassifier(mConfig, mSampleDataCache));

        if(mConfig.runClassifier(FEATURE))
            mClassifiers.add(new FeatureClassifier(mConfig, mSampleDataCache));

        if(mConfig.runClassifier(SAMPLE_TRAIT))
            mClassifiers.add(new SampleTraitClassifier(mConfig, mSampleDataCache));

        if(mConfig.runClassifier(SV))
            mClassifiers.add(new SvClassifier(mConfig, mSampleDataCache));

        if(mConfig.runClassifier(GENE_EXP))
            mClassifiers.add(new GeneExpressionClassifier(mConfig, mSampleDataCache, cmd));

        if(mConfig.runClassifier(ALT_SJ))
            mClassifiers.add(new AltSjClassifier(mConfig, mSampleDataCache));

        mSampleDataWriter = null;
        mSampleSimilarityWriter = null;
    }

    private void loadSampleData(final CommandLine cmd)
    {
        mSampleDataCache.loadReferenceSampleData(mConfig.RefSampleDataFile);
        mSampleDataCache.loadSampleData(cmd.getOptionValue(SPECIFIC_SAMPLE_DATA), mConfig.SampleDataFile);
        mSampleDataCache.markRefSamples();
    }

    public void run()
    {
        if(!mConfig.isValid() || !mSampleDataCache.isValid())
        {
            CUP_LOGGER.error("invalid config");
            return;
        }

        if(mSampleDataCache.SampleIds.isEmpty())
        {
            CUP_LOGGER.error("no samples specified");
            return;
        }

        if(!allClassifiersValid())
            return;

        if(mSampleDataCache.isMultiSample())
        {
            CUP_LOGGER.info("loaded {} samples, {} ref samples and {} ref cancer types",
                    mSampleDataCache.SampleIds.size(), mSampleDataCache.RefSampleCancerTypeMap.size(), mSampleDataCache.RefCancerSampleData.size());
        }

        initialiseOutputFiles();

        if(mSampleDataCache.SpecificSample != null)
        {
            final SampleData specificSample = mSampleDataCache.SpecificSample;

            CUP_LOGGER.info("sample({}) running CUP analysis", specificSample.Id);
            processSample(specificSample);
        }
        else
        {
            int sampleCount = 0;
            for(SampleData sample : mSampleDataCache.SampleDataList)
            {
                CUP_LOGGER.debug("sample({}) running CUP analysis", sample.Id);

                processSample(sample);

                if(!allClassifiersValid())
                    break;

                ++sampleCount;

                if((sampleCount % 100) == 0)
                {
                    CUP_LOGGER.info("processed {} samples", sampleCount);
                }
            }
        }

        closeBufferedWriter(mSampleDataWriter);
        closeBufferedWriter(mSampleSimilarityWriter);

        CUP_LOGGER.info("CUP analysis complete {}", !allClassifiersValid() ? "with errors" : "");
    }

    private boolean allClassifiersValid()
    {
        boolean allInvalid = true;
        for(CuppaClassifier classifier : mClassifiers)
        {
            if(!classifier.isValid())
            {
                allInvalid = false;
                CUP_LOGGER.error("invalid classifier({})", classifier.categoryType());
            }
        }

        return allInvalid;
    }

    private void processSample(final SampleData sample)
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

        if(mConfig.runClassifier(GENE_EXP) && (mConfig.IncludedCategories.isEmpty() || mConfig.IncludedCategories.size() > 1))
        {
            final List<SampleResult> dnaResults = allResults.stream()
                    .filter(x -> x.Category == CLASSIFIER && isDna(ClassifierType.valueOf(x.DataType)))
                    .collect(Collectors.toList());

            SampleResult dnaScoreResult = calcClassifierScoreResult(sample, dnaResults, "DNA_COMBINED");

            final List<SampleResult> rnaResults = allResults.stream()
                    .filter(x -> x.Category == CLASSIFIER && isRna(ClassifierType.valueOf(x.DataType)))
                    .collect(Collectors.toList());

            SampleResult rnaScoreResult = calcClassifierScoreResult(sample, rnaResults, "RNA_COMBINED");

            if(dnaScoreResult != null)
                allResults.add(dnaScoreResult);

            if(rnaScoreResult != null)
                allResults.add(rnaScoreResult);
        }

        SampleResult classifierScoreResult = calcClassifierScoreResult(sample, allResults, COMBINED.toString());

        if(classifierScoreResult != null)
            allResults.add(classifierScoreResult);

        writeSampleData(sample, allResults);
        writeSampleSimilarities(sample, similarities);
    }

    private void initialiseOutputFiles()
    {
        try
        {
            final String sampleDataFilename = mSampleDataCache.isSingleSample() ?
                    mConfig.OutputDir + mSampleDataCache.SpecificSample.Id + ".cup.data.csv"
                    : mConfig.formOutputFilename("SAMPLE_DATA");

            mSampleDataWriter = createBufferedWriter(sampleDataFilename, false);

            mSampleDataWriter.write("SampleId,Category,ResultType,DataType,Value,RefCancerType,RefValue");

            mSampleDataWriter.newLine();

            if(mConfig.WriteSimilarities)
            {
                final String sampleSimilarityFilename = mSampleDataCache.isSingleSample() ?
                        mConfig.OutputDir + mSampleDataCache.SpecificSample.Id + ".cup.similarities.csv"
                        : mConfig.formOutputFilename("SAMPLE_SIMILARITIES");

                mSampleSimilarityWriter = createBufferedWriter(sampleSimilarityFilename, false);

                mSampleSimilarityWriter.write("SampleId,MatchType,Score,MatchSampleId,MatchCancerType");
                mSampleSimilarityWriter.write(",MatchPrimaryType,MatchPrimarySubtype,MatchLocation,MatchSubLocation");
                mSampleSimilarityWriter.newLine();
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write CUPPA output: {}", e.toString());
        }
    }

    private void writeSampleData(final SampleData sampleData, final List<SampleResult> results)
    {
        if(results.isEmpty() || mSampleDataWriter == null)
            return;

        try
        {
            for(SampleResult result : results)
            {
                if(mConfig.WriteClassifiersOnly && result.Category != CLASSIFIER)
                    continue;

                final String sampleStr = String.format("%s,%s,%s,%s,%s",
                        sampleData.Id, result.Category, result.ResultType, result.DataType, result.Value.toString());

                for(Map.Entry<String,Double> cancerValues : result.CancerTypeValues.entrySet())
                {
                    mSampleDataWriter.write(String.format("%s,%s,%.3g",
                            sampleStr, cancerValues.getKey(), cancerValues.getValue()));
                    mSampleDataWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to write sample data: {}", e.toString());
        }
    }

    private void writeSampleSimilarities(final SampleData sampleData, final List<SampleSimilarity> similarities)
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

                mSampleSimilarityWriter.write(String.format("%s,%s,%.3f,%s",
                        sampleData.Id, similarity.MatchType, similarity.Score, similarity.MatchedSampleId));

                if(matchedSample != null)
                {
                    mSampleSimilarityWriter.write(String.format(",%s,%s,%s,%s,%s",
                            matchedSample.CancerType, matchedSample.PrimaryType, matchedSample.PrimarySubtype,
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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        CuppaConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if (cmd.hasOption(LOG_DEBUG))
        {
            Configurator.setRootLevel(Level.DEBUG);
        }

        CupAnalyser cupAnalyser = new CupAnalyser(cmd);
        cupAnalyser.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}

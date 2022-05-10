package com.hartwig.hmftools.cup;

import static java.lang.Math.min;

import static com.hartwig.hmftools.cup.CuppaConfig.LOG_DEBUG;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.SPECIFIC_SAMPLE_DATA;
import static com.hartwig.hmftools.cup.common.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.cup.common.CategoryType.FEATURE;
import static com.hartwig.hmftools.cup.common.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.cup.common.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.cup.common.CategoryType.SNV;
import static com.hartwig.hmftools.cup.common.CategoryType.SV;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.cup.common.CategoryType;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.feature.FeatureClassifier;
import com.hartwig.hmftools.cup.rna.AltSjClassifier;
import com.hartwig.hmftools.cup.rna.GeneExpressionClassifier;
import com.hartwig.hmftools.cup.traits.SampleTraitClassifier;
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
    private final ResultsWriter mResultsWriter;

    public CupAnalyser(final CommandLine cmd)
    {
        mConfig = new CuppaConfig(cmd);

        mSampleDataCache = new SampleDataCache();

        mResultsWriter = new ResultsWriter(mConfig, mSampleDataCache);

        mClassifiers = Lists.newArrayList();

        loadSampleData(cmd);

        if(!mConfig.isValid() || !mSampleDataCache.isValid())
        {
            CUP_LOGGER.error("invalid config");
            System.exit(1);
        }

        if(mConfig.runClassifier(SAMPLE_TRAIT))
            mClassifiers.add(new SampleTraitClassifier(mConfig, mSampleDataCache));

        if(mConfig.runClassifier(SNV))
            mClassifiers.add(new SomaticClassifier(mConfig, mSampleDataCache, cmd));

        if(mConfig.runClassifier(FEATURE))
            mClassifiers.add(new FeatureClassifier(mConfig, mSampleDataCache, cmd));

        if(mConfig.runClassifier(SV))
            mClassifiers.add(new SvClassifier(mConfig, mSampleDataCache));

        if(mConfig.runClassifier(GENE_EXP))
            mClassifiers.add(new GeneExpressionClassifier(mConfig, mSampleDataCache, cmd));

        if(mConfig.runClassifier(ALT_SJ))
            mClassifiers.add(new AltSjClassifier(mConfig, mSampleDataCache, cmd));

        CUP_LOGGER.debug("{} classifiers loaded", mClassifiers.size());
    }

    private void loadSampleData(final CommandLine cmd)
    {
        mSampleDataCache.loadReferenceSampleData(mConfig.RefSampleDataFile);
        mSampleDataCache.loadSampleData(cmd.getOptionValue(SPECIFIC_SAMPLE_DATA), mConfig.SampleDataFile);

        if(mSampleDataCache.isMultiSample())
        {
            CUP_LOGGER.info("loaded samples ref({}) test({}) ref cancer types({})",
                    mSampleDataCache.RefSampleCancerTypeMap.size(), mSampleDataCache.SampleIds.size(),
                    mSampleDataCache.RefCancerSampleData.size());

            if(mConfig.Categories.stream().anyMatch(x -> CategoryType.isRna(x)))
            {
                CUP_LOGGER.info("RNA samples: ref({}) test({})",
                        mSampleDataCache.RefSampleDataList.stream().filter(x -> x.hasRna()).count(),
                        mSampleDataCache.SampleDataList.stream().filter(x -> x.hasRna()).count());
            }
        }
    }

    public void run()
    {
        if(mSampleDataCache.SampleIds.isEmpty())
        {
            CUP_LOGGER.error("no samples specified");
            return;
        }

        if(!allClassifiersValid(mClassifiers))
        {
            System.exit(1);
            return;
        }

        if(mSampleDataCache.SpecificSample != null)
        {
            final SampleData specificSample = mSampleDataCache.SpecificSample;

            CUP_LOGGER.info("sample({}) running CUP analysis", specificSample.Id);
            SampleTask sampleTask = new SampleTask(0, mConfig, mSampleDataCache, mClassifiers, mResultsWriter);
            sampleTask.processSample(specificSample);
        }
        else
        {
            List<SampleTask> sampleTasks = Lists.newArrayList();

            for(int i = 0; i < min(mConfig.Threads, mSampleDataCache.SampleDataList.size()); ++i)
            {
                sampleTasks.add(new SampleTask(i, mConfig, mSampleDataCache, mClassifiers, mResultsWriter));
            }

            int taskIndex = 0;
            for(SampleData sample : mSampleDataCache.SampleDataList)
            {
                sampleTasks.get(taskIndex).getSamples().add(sample);
                ++taskIndex;

                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;
            }

            List<Callable> callableTasks = sampleTasks.stream().collect(Collectors.toList());

            TaskExecutor.executeTasks(callableTasks, mConfig.Threads);
        }

        mResultsWriter.close();

        mClassifiers.forEach(x -> x.close());

        if(!allClassifiersValid(mClassifiers))
        {
            CUP_LOGGER.info("CUP exiting with errors");
            System.exit(1);
        }

        CUP_LOGGER.info("CUP analysis complete");
    }

    public synchronized static boolean allClassifiersValid(final List<CuppaClassifier> classifiers)
    {
        boolean allInvalid = true;
        for(CuppaClassifier classifier : classifiers)
        {
            if(!classifier.isValid())
            {
                allInvalid = false;
                CUP_LOGGER.error("invalid classifier({})", classifier.categoryType());
            }
        }

        return allInvalid;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        Options options = new Options();
        CuppaConfig.addCmdLineArgs(options);

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
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

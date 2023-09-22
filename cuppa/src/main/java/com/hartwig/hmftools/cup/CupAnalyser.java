package com.hartwig.hmftools.cup;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.SAMPLE_RNA_LENGTH;
import static com.hartwig.hmftools.common.cuppa.CategoryType.ALT_SJ;
import static com.hartwig.hmftools.common.cuppa.CategoryType.FEATURE;
import static com.hartwig.hmftools.common.cuppa.CategoryType.GENE_EXP;
import static com.hartwig.hmftools.common.cuppa.CategoryType.SAMPLE_TRAIT;
import static com.hartwig.hmftools.common.cuppa.CategoryType.SNV;
import static com.hartwig.hmftools.common.cuppa.CategoryType.SV;
import static com.hartwig.hmftools.cup.common.CupConstants.APP_NAME;
import static com.hartwig.hmftools.cup.common.CupConstants.DEFAULT_RNA_LENGTH;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.r.RExecutor;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.cup.common.CuppaClassifier;
import com.hartwig.hmftools.cup.common.SampleData;
import com.hartwig.hmftools.cup.common.SampleDataCache;
import com.hartwig.hmftools.cup.feature.FeatureClassifier;
import com.hartwig.hmftools.cup.rna.AltSjClassifier;
import com.hartwig.hmftools.cup.rna.GeneExpressionClassifier;
import com.hartwig.hmftools.cup.traits.SampleTraitClassifier;
import com.hartwig.hmftools.cup.somatics.SomaticClassifier;
import com.hartwig.hmftools.cup.svs.SvClassifier;

import org.jetbrains.annotations.NotNull;

public class CupAnalyser
{
    private final CuppaConfig mConfig;

    private final SampleDataCache mSampleDataCache;

    private final List<CuppaClassifier> mClassifiers;
    private final ResultsWriter mResultsWriter;

    public CupAnalyser(final ConfigBuilder configBuilder)
    {
        mConfig = new CuppaConfig(configBuilder);

        mSampleDataCache = new SampleDataCache();

        mClassifiers = Lists.newArrayList();

        loadSampleData(configBuilder);

        if(!mConfig.isValid() || !mSampleDataCache.isValid())
        {
            CUP_LOGGER.error("invalid config");
            System.exit(1);
        }

        if(mConfig.runClassifier(SAMPLE_TRAIT)) // added first since other classifiers depend on ploidy & purity
            mClassifiers.add(new SampleTraitClassifier(mConfig, mSampleDataCache, configBuilder));

        if(mConfig.runClassifier(SNV))
            mClassifiers.add(new SomaticClassifier(mConfig, mSampleDataCache, configBuilder));

        if(mConfig.runClassifier(FEATURE))
            mClassifiers.add(new FeatureClassifier(mConfig, mSampleDataCache, configBuilder));

        if(mConfig.runClassifier(SV))
            mClassifiers.add(new SvClassifier(mConfig, mSampleDataCache));

        if(mConfig.runClassifier(GENE_EXP))
            mClassifiers.add(new GeneExpressionClassifier(mConfig, mSampleDataCache, configBuilder));

        if(mConfig.runClassifier(ALT_SJ))
            mClassifiers.add(new AltSjClassifier(mConfig, mSampleDataCache, configBuilder));

        CUP_LOGGER.debug("{} classifiers loaded", mClassifiers.size());

        mResultsWriter = new ResultsWriter(mConfig, mSampleDataCache);
    }

    private void loadSampleData(final ConfigBuilder configBuilder)
    {
        mSampleDataCache.loadReferenceSampleData(mConfig.RefSampleDataFile);

        String sampleId = configBuilder.getValue(SAMPLE);
        int rnaReadLength = Integer.parseInt(configBuilder.getValue(SAMPLE_RNA_LENGTH, String.valueOf(DEFAULT_RNA_LENGTH)));
        mSampleDataCache.loadSampleData(sampleId, rnaReadLength, mConfig.SampleDataFile);

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
            System.exit(1);
        }

        for(CuppaClassifier classifier : mClassifiers)
        {
            if(!classifier.loadData())
            {
                CUP_LOGGER.error("classifier({}) failed data loading", classifier.categoryType());
                System.exit(1);
            }
        }

        if(mSampleDataCache.SpecificSample != null)
        {
            final SampleData specificSample = mSampleDataCache.SpecificSample;

            CUP_LOGGER.info("sample({}) running CUP analysis", specificSample.Id);
            SampleTask sampleTask = new SampleTask(0, mConfig, mSampleDataCache, mClassifiers, mResultsWriter);

            if(!sampleTask.processSample(specificSample))
                System.exit(1);
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

            if(!TaskExecutor.executeTasks(callableTasks, mConfig.Threads))
                System.exit(1);
        }

        mResultsWriter.close();

        mClassifiers.forEach(x -> x.close());

        if(mConfig.CreatePdf && mSampleDataCache.SpecificSample != null)
        {
            try
            {
                RExecutor.executeFromClasspath(
                        "r/CupGenerateReport_pipeline.R", mSampleDataCache.SpecificSample.Id, mConfig.OutputDir);
            }
            catch(Exception e)
            {
                CUP_LOGGER.error("failed to generate PDF report with R script: {}", e.toString());
                System.exit(1);
            }
        }

        CUP_LOGGER.info("CUP analysis complete");
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        CuppaConfig.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        CupAnalyser cupAnalyser = new CupAnalyser(configBuilder);
        cupAnalyser.run();
    }
}

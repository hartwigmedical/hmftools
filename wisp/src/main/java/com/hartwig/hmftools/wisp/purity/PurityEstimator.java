package com.hartwig.hmftools.wisp.purity;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.wisp.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.purity.WriteType.plotCopyNumber;
import static com.hartwig.hmftools.wisp.purity.WriteType.plotSomatics;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.FittedPurityScore;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurity;
import com.hartwig.hmftools.common.purple.ImmutableFittedPurityScore;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.RunMode;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.wisp.purity.loh.AmberLohCalcs;
import com.hartwig.hmftools.wisp.purity.loh.AmberLohResult;
import com.hartwig.hmftools.wisp.purity.cn.CnPurityResult;
import com.hartwig.hmftools.wisp.purity.cn.CopyNumberProfile;
import com.hartwig.hmftools.wisp.purity.variant.SomaticPurityResult;
import com.hartwig.hmftools.wisp.purity.variant.SomaticVariants;

public class PurityEstimator
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    public PurityEstimator(final ConfigBuilder configBuilder)
    {
        mConfig = new PurityConfig(configBuilder);

        if(mConfig.Samples.isEmpty())
            System.exit(1);

        mResultsWriter = new ResultsWriter(mConfig);
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();

        List<PurityTask> purityCalcTasks = Lists.newArrayList();
        List<PlotTask> plotTasks = Lists.newArrayList();

        boolean requirePlots = plotSomatics(mConfig.WriteTypes) || plotCopyNumber(mConfig.WriteTypes);

        if(mConfig.Threads > 1)
        {
            for(int i = 0; i < min(mConfig.Samples.size(), mConfig.Threads); ++i)
            {
                purityCalcTasks.add(new PurityTask());

                if(requirePlots)
                    plotTasks.add(new PlotTask());
            }

            int taskIndex = 0;
            for(SampleData sample : mConfig.Samples)
            {
                if(taskIndex >= purityCalcTasks.size())
                    taskIndex = 0;

                purityCalcTasks.get(taskIndex).Samples.add(sample);

                if(requirePlots)
                    plotTasks.get(taskIndex).Samples.add(sample);

                ++taskIndex;
            }

            CT_LOGGER.debug("splitting {} patients across {} threads", mConfig.Samples.size(), purityCalcTasks.size());

            List<Callable> callableList = purityCalcTasks.stream().collect(Collectors.toList());
            if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
            {
                System.exit(1);
            }
        }
        else
        {
            PurityTask sampleTask = new PurityTask();
            sampleTask.Samples.addAll(mConfig.Samples);
            sampleTask.call();

            if(requirePlots)
            {
                PlotTask plotTask = new PlotTask();
                plotTask.Samples.addAll(mConfig.Samples);
                plotTasks.add(plotTask);
            }
        }

        mResultsWriter.close();

        if(requirePlots)
        {
            CT_LOGGER.debug("generating plots");

            if(plotTasks.size() > 1)
            {
                final List<Callable> callableList = plotTasks.stream().collect(Collectors.toList());
                if(!TaskExecutor.executeTasks(callableList, mConfig.Threads))
                {
                    System.exit(1);
                }
            }
            else
            {
                plotTasks.get(0).call();
            }
        }

        CT_LOGGER.info("Wisp purity estimator complete{}",
                mConfig.Samples.size() > 1 ? format(", mins(%s)", runTimeMinsStr(startTimeMs)) : "");
    }

    private class PurityTask implements Callable
    {
        public final List<SampleData> Samples;

        public PurityTask()
        {
            Samples = Lists.newArrayList();
        }

        @Override
        public Long call()
        {
            for(SampleData sample : Samples)
            {
                try
                {
                    processSample(sample);
                }
                catch(Throwable t)
                {
                    CT_LOGGER.error("failed processing patient({}): {}", sample.PatientId, t.toString());
                    t.printStackTrace();
                    System.exit(1);
                }
            }

            return (long)0;
        }

        private void processSample(final SampleData sample)
        {
            CT_LOGGER.info("processing sample: {}", sample);

            PurityContext purityContext = loadPurplePurity(sample);

            SomaticVariants somaticVariants = null;

            if(mConfig.PurityMethods.contains(PurityMethod.SOMATIC_VARIANT))
            {
                somaticVariants = new SomaticVariants(mConfig, mResultsWriter, sample);
                if(!somaticVariants.loadVariants())
                {
                    if(!mConfig.AllowMissingSamples)
                        System.exit(1);
                }
            }

            CopyNumberProfile copyNumberProfile = null;
            if(!sample.IsPanel && mConfig.PurityMethods.contains(PurityMethod.COPY_NUMBER))
            {
                copyNumberProfile = new CopyNumberProfile(mConfig, mResultsWriter, sample);

                if(!copyNumberProfile.hasValidData())
                {
                    if(!mConfig.AllowMissingSamples)
                        System.exit(1);

                    CT_LOGGER.warn("sample({}) has missing Cobalt data", sample);
                }
            }

            AmberLohCalcs amberLohCalcs = null;
            if(!sample.IsPanel && mConfig.PurityMethods.contains(PurityMethod.AMBER_LOH))
            {
                amberLohCalcs = new AmberLohCalcs(mConfig, mResultsWriter, sample);

                if(!amberLohCalcs.hasValidData())
                {
                    if(!mConfig.AllowMissingSamples)
                        System.exit(1);

                    CT_LOGGER.warn("sample({}) has missing Amber data", sample);
                }
            }

            for(String sampleId : sample.SampleIds)
            {
                CnPurityResult cnPurityResult = copyNumberProfile != null ?
                        copyNumberProfile.processSample(sampleId, purityContext) : CnPurityResult.INVALID_RESULT;

                SomaticPurityResult somaticPurityResult = somaticVariants != null ?
                        somaticVariants.processSample(sampleId, purityContext) : SomaticPurityResult.INVALID_RESULT;

                AmberLohResult lohResult = amberLohCalcs != null ? amberLohCalcs.processSample(sampleId) : AmberLohResult.INVALID_RESULT;

                if(cnPurityResult == null || somaticPurityResult == null || lohResult == null)
                {
                    if(!mConfig.AllowMissingSamples)
                    {
                        System.exit(1);
                    }
                    else
                    {
                        CT_LOGGER.warn("sample({}) has missing data", sample);
                    }
                }

                mResultsWriter.writeSampleSummary(sample, sampleId, purityContext, cnPurityResult, somaticPurityResult, lohResult);
            }
        }

        private PurityContext loadPurplePurity(final SampleData sample)
        {
            if(!mConfig.hasSyntheticTumor())
            {
                try
                {
                    return PurityContextFile.read(mConfig.getPurpleDir(sample.TumorId), sample.TumorId);
                }
                catch(Exception e)
                {
                    CT_LOGGER.error("failed to load Purple purity: {}", e.toString());
                    System.exit(1);
                }
            }
            else
            {
                double purity = 1;
                double ploidy = 2;

                FittedPurity fittedPurity = ImmutableFittedPurity.builder()
                        .purity(purity).normFactor(1).ploidy(ploidy).score(0D).diploidProportion(0D).somaticPenalty(0D)
                        .build();

                PurpleQC purpleQC = ImmutablePurpleQC.builder()
                        .method(FittedPurityMethod.NORMAL).amberMeanDepth(0).copyNumberSegments(1).unsupportedCopyNumberSegments(0)
                        .deletedGenes(0).purity(purity).contamination(0D).cobaltGender(Gender.FEMALE).amberGender(Gender.FEMALE)
                        .lohPercent(0)
                        .build();

                FittedPurityScore fittedPurityScore = ImmutableFittedPurityScore.builder()
                        .minPurity(0D).maxPurity(0D).minPloidy(0D).maxPloidy(0D).minDiploidProportion(0D).maxDiploidProportion(0D)
                        .build();

                PurityContext genericContext = ImmutablePurityContext.builder()
                        .gender(Gender.FEMALE).runMode(RunMode.TUMOR_GERMLINE).targeted(false).bestFit(fittedPurity)
                        .method(FittedPurityMethod.NORMAL).score(fittedPurityScore).qc(purpleQC).polyClonalProportion(0D)
                        .wholeGenomeDuplication(false).microsatelliteIndelsPerMb(0D).tumorMutationalBurdenPerMb(0D)
                        .tumorMutationalLoad(0).svTumorMutationalBurden(0).microsatelliteStatus(MicrosatelliteStatus.UNKNOWN)
                        .tumorMutationalLoadStatus(TumorMutationalStatus.UNKNOWN).tumorMutationalBurdenStatus(TumorMutationalStatus.UNKNOWN)
                        .build();

                return genericContext;
            }

            return null;
        }
    }

    private class PlotTask implements Callable
    {
        public final List<SampleData> Samples;

        public PlotTask()
        {
            Samples = Lists.newArrayList();
        }

        @Override
        public Long call()
        {
            for(SampleData sample : Samples)
            {
                createSamplePlots(sample);
            }

            return (long)0;
        }

        private void createSamplePlots(final SampleData sample)
        {
            // CT_LOGGER.info("processing sample: {}", sample);

            for(String sampleId : sample.SampleIds)
            {
                if(!sample.IsPanel && plotCopyNumber(mConfig.WriteTypes))
                {
                    if(!CopyNumberProfile.plotCopyNumberGcRatioFit(sample.PatientId, sampleId, mConfig))
                        return;
                }

                if(plotSomatics(mConfig.WriteTypes))
                {
                    if(!SomaticVariants.plotSomaticVafs(sample.PatientId, sampleId, mConfig))
                    {
                        CT_LOGGER.error("patient({}) sample({}) somatic plot failed", sample.PatientId, sampleId);
                    }
                }
            }
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        PurityConfig.addConfig(configBuilder);

        ConfigUtils.addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        PurityEstimator purityEstimator = new PurityEstimator(configBuilder);
        purityEstimator.run();
    }
}

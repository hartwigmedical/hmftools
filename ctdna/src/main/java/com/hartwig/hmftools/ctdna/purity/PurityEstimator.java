package com.hartwig.hmftools.ctdna.purity;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.ctdna.purity.cn.CnPurityResult;
import com.hartwig.hmftools.ctdna.purity.cn.CopyNumberProfile;
import com.hartwig.hmftools.ctdna.purity.variant.SomaticVariantResult;
import com.hartwig.hmftools.ctdna.purity.variant.SomaticVariants;

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
        List<PurityTask> sampleTasks = Lists.newArrayList();

        if(mConfig.Threads > 1)
        {
            for(int i = 0; i < min(mConfig.Samples.size(), mConfig.Threads); ++i)
            {
                sampleTasks.add(new PurityTask());
            }

            int taskIndex = 0;
            for(SampleData sample : mConfig.Samples)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).Samples.add(sample);

                ++taskIndex;
            }

            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            PurityTask sampleTask = new PurityTask();
            sampleTask.Samples.addAll(mConfig.Samples);
            sampleTask.call();
        }

        mResultsWriter.close();

        if(mConfig.PlotCnFit && mConfig.WriteCnRatios)
        {
            boolean hasError = false;
            for(SampleData sample : mConfig.Samples)
            {
                for(String sampleId : sample.CtDnaSamples)
                {
                    if(!CopyNumberProfile.plotCopyNumberGcRatioFit(sample.PatientId, sampleId, mConfig))
                    {
                        hasError = true;
                        break;
                    }
                }

                if(hasError)
                    break;
            }
        }

        CT_LOGGER.info("CtDNA purity estimator complete");
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
                processSample(sample);
            }

            return (long)0;
        }

        private void processSample(final SampleData sample)
        {
            CT_LOGGER.info("processing sample: {}", sample);

            PurityContext purityContext = null;

            try
            {
                purityContext = PurityContextFile.read(mConfig.PurpleDir, sample.TumorId);
            }
            catch(Exception e)
            {
                CT_LOGGER.error("failed to load Purple purity: {}", e.toString());
                System.exit(1);
            }

            SomaticVariants somaticVariants = null;

            if(mConfig.PurityMethods.contains(PurityMethod.SOMATIC))
            {
                somaticVariants = new SomaticVariants(mConfig, mResultsWriter, sample);
                if(!somaticVariants.loadVariants())
                    System.exit(1);
            }

            CopyNumberProfile copyNumberProfile = null;
            if(mConfig.PurityMethods.contains(PurityMethod.COPY_NUMBER))
            {
                copyNumberProfile = new CopyNumberProfile(mConfig, mResultsWriter, sample);
            }

            for(String ctDnaSample : sample.CtDnaSamples)
            {
                CnPurityResult cnPurityResult = copyNumberProfile != null ?
                        copyNumberProfile.processSample(ctDnaSample, purityContext) : CnPurityResult.INVALID_RESULT;

                SomaticVariantResult somaticVariantResult = somaticVariants != null ?
                        somaticVariants.processSample(ctDnaSample, purityContext) : SomaticVariantResult.INVALID_RESULT;

                mResultsWriter.writeSampleSummary(sample.PatientId, ctDnaSample, purityContext, cnPurityResult, somaticVariantResult);
            }
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PurityConfig.addConfig(configBuilder);

        ConfigUtils.addLoggingOptions(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);

        PurityEstimator purityEstimator = new PurityEstimator(configBuilder);
        purityEstimator.run();
    }
}

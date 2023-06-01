package com.hartwig.hmftools.ctdna.purity;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;

import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class PurityEstimator
{
    private final PurityConfig mConfig;
    private final ResultsWriter mResultsWriter;

    public PurityEstimator(final CommandLine cmd)
    {
        mConfig = new PurityConfig(cmd);
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

    public static void main(final String[] args) throws ParseException
    {
        final Options options = new Options();
        PurityConfig.addCommandLineOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PurityEstimator purityEstimator = new PurityEstimator(cmd);
        purityEstimator.run();

        CT_LOGGER.info("CtDNA purity estimator complete");
    }

    private static CommandLine createCommandLine(final String[] args, final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}

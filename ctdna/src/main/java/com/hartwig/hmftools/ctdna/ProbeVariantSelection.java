package com.hartwig.hmftools.ctdna;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;
import static com.hartwig.hmftools.ctdna.PvConfig.createCmdLineOptions;
import static com.hartwig.hmftools.ctdna.VariantUtils.calcGcPercent;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class ProbeVariantSelection
{
    private final PvConfig mConfig;
    private final List<Variant> mCommonVariants;
    private final RefGenomeInterface mRefGenome;
    private final BufferedWriter mWriter;

    public ProbeVariantSelection(final CommandLine cmd)
    {
        mConfig = new PvConfig(cmd);
        mRefGenome = loadRefGenome(mConfig.RefGenomeFile);

        mCommonVariants = Lists.newArrayList();

        mWriter = initialiseWriter();
    }

    public void run()
    {
        if(!mConfig.isValid() || mWriter == null || mRefGenome == null)
            System.exit(1);

        if(mConfig.isMultiSample())
            PV_LOGGER.info("running probe variant selection for {} samples", mConfig.SampleIds.size());
        else
            PV_LOGGER.info("sample({}) running probe variant selection", mConfig.sample());

        if(mConfig.ReferenceVariantsFile != null)
        {
            List<Variant> referenceVariants = ReferenceMutation.loadKnownMutations(mConfig.ReferenceVariantsFile);

            if(referenceVariants == null)
                System.exit(1);

            mCommonVariants.addAll(referenceVariants);
        }

        List<SampleTask> sampleTasks = Lists.newArrayList();

        if(mConfig.Threads > 1)
        {
            for(int i = 0; i < min(mConfig.SampleIds.size(), mConfig.Threads); ++i)
            {
                sampleTasks.add(new SampleTask(i));
            }

            int taskIndex = 0;
            for(String sampleId : mConfig.SampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

                ++taskIndex;
            }

            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            SampleTask sampleTask = new SampleTask(0);
            sampleTask.getSampleIds().addAll(mConfig.SampleIds);
            sampleTasks.add(sampleTask);
            sampleTask.call();
        }

        closeBufferedWriter(mWriter);

        PV_LOGGER.info("Probe variation selection complete");
    }

    private class SampleTask implements Callable
    {
        private final int mTaskId;
        private final List<String> mSampleIds;

        public SampleTask(int taskId)
        {
            mTaskId = taskId;
            mSampleIds = Lists.newArrayList();
        }

        public List<String> getSampleIds() { return mSampleIds; }

        @Override
        public Long call()
        {
            for(int i = 0; i < mSampleIds.size(); ++i)
            {
                String sampleId = mSampleIds.get(i);

                processSample(sampleId);
                if(i > 0 && (i % 10) == 0)
                {
                    PV_LOGGER.info("{}: processed {} samples", mTaskId, i);
                }
            }

            if(mConfig.Threads > 1)
            {
                PV_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());
            }

            return (long)0;
        }

        private void processSample(final String sampleId)
        {
            List<Variant> variants = Lists.newArrayList();
            variants.addAll(mCommonVariants);
            variants.addAll(PointMutation.loadSomatics(sampleId, mConfig));
            variants.addAll(GermlineMutation.loadGermlineMutations(sampleId, mConfig));
            variants.addAll(StructuralVariant.loadStructuralVariants(sampleId, mConfig));
            variants.addAll(GermlineSv.loadGermlineStructuralVariants(sampleId, mConfig));

            variants.forEach(x -> x.generateSequences(mRefGenome, mConfig));

            List<Variant> selectedVariants = VariantSelection.selectVariants(variants, mConfig);

            writeSelectedVariants(sampleId, selectedVariants);
        }
    }

    private BufferedWriter initialiseWriter()
    {
        if(mConfig.SampleIds.isEmpty())
            return null;

        try
        {
            String filename = mConfig.OutputDir;

            if(mConfig.isMultiSample())
                filename += "cohort_probe_variants.csv";
            else
                filename += mConfig.sample() + ".probe_variants.csv";

            BufferedWriter writer = createBufferedWriter(filename, false);

            if(mConfig.isMultiSample())
                writer.write("SampleId,");

            writer.write("Category,Variant,CopyNumber,Vaf,TumorFrags,PhasedVariants,Gene");
            writer.write(",Type,Sequence,GcPercent");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            PV_LOGGER.error(" failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private synchronized void writeSelectedVariants(final String sampleId, final List<Variant> selectedVariants)
    {
        if(mWriter == null)
            return;

        try
        {
            for(Variant variant : selectedVariants)
            {
                String variantInfo = mConfig.isMultiSample() ? format("%s,", sampleId) : "";

                variantInfo += format("%s,%s,%.2f,%.2f,%d,%s,%s",
                        variant.categoryType(), variant.description(), variant.copyNumber(), variant.vaf(),
                        variant.tumorFragments(), variant.hasPhaseVariants(), variant.gene());

                mWriter.write(format("%s,%s,%s,%.2f", variantInfo, "ALT", variant.sequence(), variant.gc()));
                mWriter.newLine();

                for(String refSequence : variant.refSequences())
                {
                    mWriter.write(format("%s,%s,%s,%.2f", variantInfo, "REF", refSequence, calcGcPercent(refSequence)));
                    mWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error(" failed to write variants: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        final VersionInfo version = new VersionInfo("ctdna.version");
        PV_LOGGER.info("ProbeVariantSelection version: {}", version.version());

        final Options options = createCmdLineOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            ProbeVariantSelection application = new ProbeVariantSelection(cmd);
            application.run();
        }
        catch(ParseException e)
        {
            PV_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PrimerVariantSelection", options);
            System.exit(1);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}

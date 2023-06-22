package com.hartwig.hmftools.ctdna.probe;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.ctdna.common.CommonUtils.calcGcPercent;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.jetbrains.annotations.NotNull;

public class ProbeFinder
{
    private final PvConfig mConfig;
    private final List<Variant> mCommonVariants;
    private final RefGenomeInterface mRefGenome;
    private final BufferedWriter mWriter;

    public ProbeFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new PvConfig(configBuilder);
        mRefGenome = loadRefGenome(mConfig.RefGenomeFile);

        mCommonVariants = Lists.newArrayList();

        mWriter = initialiseWriter();
    }

    public void run()
    {
        if(!mConfig.isValid() || mWriter == null || mRefGenome == null)
            System.exit(1);

        if(mConfig.isMultiSample())
            CT_LOGGER.info("running probe variant selection for {} samples", mConfig.SampleIds.size());
        else
            CT_LOGGER.info("sample({}) running probe variant selection", mConfig.sample());

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

        CT_LOGGER.info("Probe variation selection complete");
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
                    CT_LOGGER.info("{}: processed {} samples", mTaskId, i);
                }
            }

            if(mConfig.Threads > 1)
            {
                CT_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());
            }

            return (long)0;
        }

        private void processSample(final String sampleId)
        {
            try
            {
                List<Variant> variants = Lists.newArrayList();
                variants.addAll(mCommonVariants);
                variants.addAll(PointMutation.loadSomatics(sampleId, mConfig));
                variants.addAll(GermlineMutation.loadGermlineMutations(sampleId, mConfig));
                variants.addAll(StructuralVariant.loadStructuralVariants(sampleId, mConfig));
                variants.addAll(GermlineSv.loadGermlineStructuralVariants(sampleId, mConfig));

                variants.forEach(x -> x.generateSequences(mRefGenome, mConfig));

                List<Variant> selectedVariants = VariantSelection.selectVariants(variants, mConfig);

                writeVariants(sampleId, selectedVariants);
            }
            catch(Exception e)
            {
                CT_LOGGER.error("sample data loading error: {}", e.toString());
                e.printStackTrace();
                System.exit(1);
            }
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
                filename += "cohort_probe_variants";
            else
                filename += mConfig.sample() + ".probe_variants";

            if(mConfig.OutputId != null)
            {
                filename += "." + mConfig.OutputId;
            }

            filename += TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(mConfig.isMultiSample())
                sj.add("SampleId");

            sj.add("Category").add("Status").add("Variant").add("Reported").add("CopyNumber").add("Vaf").add("TumorFrags");
            sj.add("PhasedVariants").add("Gene").add("Type").add("Sequence").add("GcPercent").add("OtherData");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            CT_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private synchronized void writeVariants(final String sampleId, final List<Variant> selectedVariants)
    {
        if(mWriter == null)
            return;

        try
        {
            for(Variant variant : selectedVariants)
            {
                if(!mConfig.WriteAll && !variant.isSelected())
                    continue;

                StringJoiner variantInfo = new StringJoiner(TSV_DELIM);

                if(mConfig.isMultiSample())
                    variantInfo.add(sampleId);

                variantInfo.add(variant.categoryType().toString());
                variantInfo.add(variant.selectionStatus().toString());
                variantInfo.add(variant.description());
                variantInfo.add(String.valueOf(variant.reported()));
                variantInfo.add(format("%.2f", variant.copyNumber()));
                variantInfo.add(format("%.2f", variant.vaf()));
                variantInfo.add(String.valueOf(variant.tumorFragments()));
                variantInfo.add(String.valueOf(variant.hasPhaseVariants()));
                variantInfo.add(variant.gene());

                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(variantInfo.toString());
                sj.add("ALT");
                sj.add(variant.sequence());
                sj.add(format("%.2f", variant.gc()));
                sj.add(variant.otherData());

                mWriter.write(sj.toString());
                mWriter.newLine();

                for(String refSequence : variant.refSequences())
                {
                    StringJoiner refSj = new StringJoiner(TSV_DELIM);
                    refSj.add(variantInfo.toString());
                    refSj.add("REF");
                    refSj.add(refSequence);
                    refSj.add(format("%.2f", calcGcPercent(refSequence)));
                    refSj.add("");
                    mWriter.write(refSj.toString());
                    mWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            CT_LOGGER.error(" failed to write variants: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        final VersionInfo version = new VersionInfo("ctdna.version");
        CT_LOGGER.info("ProbeVariantSelection version: {}", version.version());

        ConfigBuilder configBuilder = new ConfigBuilder();
        PvConfig.addConfig(configBuilder);

        addLoggingOptions(configBuilder);

        if(!configBuilder.parseCommandLine(args))
        {
            configBuilder.logInvalidDetails();
            System.exit(1);
        }

        setLogLevel(configBuilder);

        ProbeFinder application = new ProbeFinder(configBuilder);
        application.run();
    }
}

package com.hartwig.hmftools.wisp.probe;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.wisp.common.CommonUtils.CT_LOGGER;
import static com.hartwig.hmftools.wisp.common.CommonUtils.FLD_BATCH_ID;
import static com.hartwig.hmftools.wisp.common.CommonUtils.FLD_CATEGORY;
import static com.hartwig.hmftools.wisp.common.CommonUtils.FLD_VARIANT;
import static com.hartwig.hmftools.wisp.probe.ProbeConfig.NO_BATCH_ID;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ProbeFinder
{
    private final ProbeConfig mConfig;
    private final List<Variant> mCommonVariants;
    private final RefGenomeInterface mRefGenome;
    private final Map<String, BufferedWriter> mWriters;

    public ProbeFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new ProbeConfig(configBuilder);
        mRefGenome = loadRefGenome(mConfig.RefGenomeFile);

        mCommonVariants = Lists.newArrayList();

        mWriters = Maps.newHashMap();
    }

    public void run()
    {
        if(!mConfig.isValid() || mRefGenome == null)
            System.exit(1);

        // create the batch writers
        if(mConfig.isMultiSample())
        {
            for(String batchId : mConfig.BatchSampleIds.keySet())
            {
                mWriters.put(batchId, initialiseWriter(batchId));
            }
        }
        else
        {
            mWriters.put(NO_BATCH_ID, initialiseWriter(null));
        }

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

        if(mConfig.Threads > 1 || mConfig.BatchSampleIds.size() > 1)
        {
            // split samples by batches if configured as such
            int requiredSampleTasks;

            if(mConfig.BatchSampleIds.size() <= 1)
            {
                requiredSampleTasks = min(mConfig.SampleIds.size(), mConfig.Threads);
            }
            else
            {
                requiredSampleTasks = min(mConfig.BatchSampleIds.size(), mConfig.Threads);
            }

            for(int i = 0; i < requiredSampleTasks; ++i)
            {
                sampleTasks.add(new SampleTask(i));
            }

            int taskIndex = 0;
            if(mConfig.BatchSampleIds.size() <= 1)
            {
                for(String sampleId : mConfig.SampleIds)
                {
                    if(taskIndex >= sampleTasks.size())
                        taskIndex = 0;

                    sampleTasks.get(taskIndex).addSampleId(sampleId);

                    ++taskIndex;
                }
            }
            else
            {
                for(String batchId : mConfig.BatchSampleIds.keySet())
                {
                    if(taskIndex >= sampleTasks.size())
                        taskIndex = 0;

                    sampleTasks.get(taskIndex).addBatchId(batchId);

                    ++taskIndex;
                }
            }

            final List<Callable<Void>> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mConfig.Threads);
        }
        else
        {
            SampleTask sampleTask = new SampleTask(0);
            mConfig.SampleIds.forEach(x -> sampleTask.addSampleId(x));
            sampleTasks.add(sampleTask);
            sampleTask.call();
        }

        mWriters.values().forEach(x -> closeBufferedWriter(x));

        CT_LOGGER.info("Probe variation selection complete");
    }

    private class SampleTask implements Callable<Void>
    {
        private final int mTaskId;
        private final List<String> mSampleIds;
        private final List<String> mBatchIds;

        public SampleTask(int taskId)
        {
            mTaskId = taskId;
            mSampleIds = Lists.newArrayList();
            mBatchIds = Lists.newArrayList();
        }

        public void addSampleId(final String sampleId) { mSampleIds.add(sampleId); }
        public void addBatchId(final String batchId) { mBatchIds.add(batchId); }

        @Override
        public Void call()
        {
            if(!mBatchIds.isEmpty())
            {
                for(String batchId : mBatchIds)
                {
                    Map<String, List<Variant>> sampleVariantMap = Maps.newHashMap();

                    List<String> sampleIds = mConfig.BatchSampleIds.get(batchId);

                    for(String sampleId : sampleIds)
                    {
                        List<Variant> selectedVariants = selectSampleVariants(sampleId);
                        sampleVariantMap.put(sampleId, selectedVariants);
                    }

                    StringJoiner sj = new StringJoiner(",");
                    sampleVariantMap.keySet().forEach(x -> sj.add(x));
                    CT_LOGGER.debug("batching {} samples: {}", sampleVariantMap.size(), sj.toString());

                    finaliseBatchVariants(sampleVariantMap);
                    writeVariants(CategoryType.REFERENCE.toString(), batchId, mCommonVariants);

                    for(Map.Entry<String, List<Variant>> entry : sampleVariantMap.entrySet())
                    {
                        writeVariants(entry.getKey(), batchId, entry.getValue());
                    }
                }
            }
            else
            {
                for(int i = 0; i < mSampleIds.size(); ++i)
                {
                    String sampleId = mSampleIds.get(i);
                    List<Variant> selectedVariants = selectSampleVariants(sampleId);
                    writeVariants(sampleId, null, selectedVariants);

                    if(i > 0 && (i % 10) == 0)
                    {
                        CT_LOGGER.info("{}: processed {} samples", mTaskId, i);
                    }
                }
            }

            if(mConfig.Threads > 1)
            {
                CT_LOGGER.info("{}: tasks complete for {} samples", mTaskId, max(mSampleIds.size(), mBatchIds.size()));
            }

            return null;
        }

        private List<Variant> selectSampleVariants(final String sampleId)
        {
            if(!mConfig.checkSampleDirectories(sampleId))
                return Collections.emptyList();

            try
            {
                List<Variant> variants = Lists.newArrayList();
                variants.addAll(mCommonVariants);
                variants.addAll(SomaticMutation.loadSomatics(sampleId, mConfig));
                variants.addAll(GermlineMutation.loadGermlineMutations(sampleId, mConfig));
                variants.addAll(StructuralVariant.loadStructuralVariants(sampleId, mConfig));
                variants.addAll(GermlineSv.loadGermlineStructuralVariants(sampleId, mConfig));

                variants.forEach(x -> x.generateSequences(mRefGenome, mConfig));

                return VariantSelection.selectVariants(variants, mConfig);
            }
            catch(Exception e)
            {
                CT_LOGGER.error("sample data loading error: {}", e.toString());
                e.printStackTrace();

                if(mConfig.AllowMissing)
                    return Collections.emptyList();

                System.exit(1);
            }

            return null;
        }

        private void finaliseBatchVariants(final Map<String, List<Variant>> sampleVariantMap)
        {
            // remove ref variants but register their locations for proximity checking
            ProximateLocations registeredLocations = new ProximateLocations();

            mCommonVariants.forEach(x -> x.checkAndRegisterLocation(registeredLocations));

            for(List<Variant> variants : sampleVariantMap.values())
            {
                int index = 0;
                while(index < variants.size())
                {
                    Variant variant = variants.get(index);

                    if(variant.categoryType() == CategoryType.REFERENCE)
                    {
                        variants.remove(index);
                    }
                    else
                    {
                        // remove based on proximity to another
                        if(variant.checkAndRegisterLocation(registeredLocations))
                            ++index;
                        else
                            variants.remove(index);
                    }
                }
            }
        }
    }

    private BufferedWriter initialiseWriter(@Nullable final String batchId)
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
                filename += "." + mConfig.OutputId;

            if(batchId != null)
                filename += "." + batchId;

            filename += TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(mConfig.isMultiSample())
                sj.add("SampleId");

            if(mConfig.BatchSampleIds.size() > 1)
                sj.add(FLD_BATCH_ID);

            sj.add(FLD_CATEGORY).add("Status").add(FLD_VARIANT).add("Reported").add("CopyNumber").add("Vaf").add("TumorFrags");
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

    private synchronized void writeVariants(
            final String sampleId, @Nullable final String batchId, final List<Variant> selectedVariants)
    {
        BufferedWriter writer = batchId != null ? mWriters.get(batchId) : mWriters.get(NO_BATCH_ID);

        if(writer == null)
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

                if(batchId != null)
                    variantInfo.add(batchId);

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

                writer.write(sj.toString());
                writer.newLine();

                for(String refSequence : variant.refSequences())
                {
                    StringJoiner refSj = new StringJoiner(TSV_DELIM);
                    refSj.add(variantInfo.toString());
                    refSj.add("REF");
                    refSj.add(refSequence);
                    refSj.add(format("%.2f", calcGcPercent(refSequence)));
                    refSj.add("");
                    writer.write(refSj.toString());
                    writer.newLine();
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
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        ProbeConfig.addConfig(configBuilder);

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ProbeFinder application = new ProbeFinder(configBuilder);
        application.run();
    }
}

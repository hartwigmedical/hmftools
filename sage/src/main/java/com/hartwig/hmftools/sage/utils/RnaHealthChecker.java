package com.hartwig.hmftools.sage.utils;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.Integers.median;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.sage.SageCommon.APP_NAME;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

public class RnaHealthChecker
{
    private final List<String> mSampleIds;

    private final String mVcfPath;
    private final String mOutputDir;

    private final int mDepthThreshold;
    private final int mSupportThreshold;
    private final String mRnaSampleSuffix;

    private final BufferedWriter mWriter;

    private final int mThreads;

    private static final String VCF_PATH = "vcf_path";
    private static final String MIN_SUPPORT_FRAGS = "min_support";
    private static final String MIN_DEPTH_FRAGS = "min_depth";
    private static final String RNA_SAMPLE_SUFFIX = "rna_sample_suffix";

    private static final int DEFAULT_MIN_DEPTH_FRAGS = 10;
    private static final int DEFAULT_MIN_SUPPORT_FRAGS = 1;
    private static final String DEFAULT_RNA_SAMPLE_SUFFIX = "_RNA";

    public RnaHealthChecker(final ConfigBuilder configBuilder)
    {
        mSampleIds = Lists.newArrayList();

        if(configBuilder.hasValue(SAMPLE))
        {
            mSampleIds.add(configBuilder.getValue(SAMPLE));
        }
        else if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            mSampleIds.addAll(loadSampleIdsFile(configBuilder));
        }
        else
        {
            SG_LOGGER.error("missing sample or sample ID file in config");
            System.exit(1);
        }

        mVcfPath = configBuilder.getValue(VCF_PATH);
        mRnaSampleSuffix = configBuilder.getValue(RNA_SAMPLE_SUFFIX, DEFAULT_RNA_SAMPLE_SUFFIX);

        if(mVcfPath == null)
        {
            SG_LOGGER.error("missing VCF path in config");
            System.exit(1);
        }

        mOutputDir = parseOutputDir(configBuilder);

        mDepthThreshold = configBuilder.getInteger(MIN_DEPTH_FRAGS);
        mSupportThreshold = configBuilder.getInteger(MIN_SUPPORT_FRAGS);
        mThreads = parseThreads(configBuilder);

        mWriter = initialiseWriter();
    }

    private static final int LOG_COUNT = 100000;

    public void run()
    {
        if(mSampleIds.size() == 1)
            SG_LOGGER.info("RNA Health-Checker for sample({})", mSampleIds.get(0));
        else
            SG_LOGGER.info("RNA Health-Checker for {} samples", mSampleIds.size());

        List<SampleTask> sampleTasks = com.google.common.collect.Lists.newArrayList();

        if(mThreads > 1)
        {
            for(int i = 0; i < min(mSampleIds.size(), mThreads); ++i)
            {
                sampleTasks.add(new SampleTask(i));
            }

            int taskIndex = 0;
            for(String sampleId : mSampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

                ++taskIndex;
            }

            final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mThreads);
        }
        else
        {
            SampleTask sampleTask = new SampleTask(0);
            sampleTask.getSampleIds().addAll(mSampleIds);
            sampleTasks.add(sampleTask);
            sampleTask.call();
        }

        closeBufferedWriter(mWriter);

        SG_LOGGER.info("RNA Health-Checker complete");
    }

    private class SampleTask implements Callable
    {
        private final int mTaskId;
        private final List<String> mSampleIds;

        public SampleTask(int taskId)
        {
            mTaskId = taskId;
            mSampleIds = com.google.common.collect.Lists.newArrayList();
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
                    SG_LOGGER.info("{}: processed {} samples", mTaskId, i);
                }
            }

            if(mThreads > 1)
            {
                SG_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());
            }

            return (long)0;
        }

        private void processSample(final String sampleId)
        {
            String vcfFile = mVcfPath.replaceAll("\\*", sampleId);

            if(!Files.exists(Paths.get(vcfFile)))
            {
                SG_LOGGER.warn("sample({}) missing VCF({})", sampleId, vcfFile);
                return;
            }

            SampleVariantSummary summaryData = new SampleVariantSummary();

            String rnaSampleId = sampleId + mRnaSampleSuffix;

            try
            {
                VcfFileReader vcfFileReader = new VcfFileReader(vcfFile);

                if(!vcfFileReader.fileValid())
                    return;

                VCFHeader vcfHeader = vcfFileReader.vcfHeader();

                if(vcfHeader.getGenotypeSamples().stream().noneMatch(x -> x.equals(sampleId))
                || vcfHeader.getGenotypeSamples().stream().noneMatch(x -> x.equals(rnaSampleId)))
                {
                    SG_LOGGER.warn("sample({}) VCF({}) missing genotype sample name", sampleId, vcfFile);
                    return;
                }

                for(VariantContext variantContext : vcfFileReader.iterator())
                {
                    if(variantContext.isFiltered())
                        continue;

                    if(VariantType.type(variantContext) != VariantType.SNP)
                        continue;

                    ++summaryData.VariantCount;

                    VariantImpact variantImpact = VariantImpactSerialiser.fromVariantContext(variantContext);

                    if(variantImpact == null || variantImpact.GeneName.isEmpty())
                        continue;

                    ++summaryData.GeneVariantCount;

                    Genotype genotype = variantContext.getGenotype(rnaSampleId);

                    summaryData.VariantFragments.add(new VariantSupport(genotype.getDP(), genotype.getAD()[1]));
                }
            }
            catch(Exception e)
            {
                SG_LOGGER.error("sample({}) error parsing VCF: {}", sampleId, e.toString());
            }

            SG_LOGGER.info("sample({}) processed {} variants from VCF({})", sampleId, summaryData.VariantCount, vcfFile);

            writeSampleSummary(sampleId, summaryData);
        }
    }

    private class VariantSupport
    {
        public final int Depth;
        public final int AlleleFragments;

        public VariantSupport(final int depth, final int alleleFragments)
        {
            Depth = depth;
            AlleleFragments = alleleFragments;
        }

        public String toString() { return format("depth(%d) alleleFrags(%d)", Depth, AlleleFragments); }
    }

    private class SampleVariantSummary
    {
        public int VariantCount;
        public int GeneVariantCount;
        public final List<VariantSupport> VariantFragments;

        public SampleVariantSummary()
        {
            VariantCount = 0;
            GeneVariantCount = 0;
            VariantFragments = Lists.newArrayList();
        }

        public double medianDepth()
        {
            List<Integer> depth = Lists.newArrayList();
            VariantFragments.stream().filter(x -> x.Depth > 0).forEach(x -> depth.add(x.Depth));
            return median(depth);
        }

        public double medianSupport()
        {
            List<Integer> support = Lists.newArrayList();
            VariantFragments.stream().filter(x -> x.Depth > 0).forEach(x -> support.add(x.AlleleFragments));
            return median(support);
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mOutputDir;

            if(mSampleIds.size() == 1)
                fileName += mSampleIds.get(0);
            else
                fileName += "cohort";

            fileName += ".rna_health_check" + TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(fileName, false);

            if(mSampleIds.size() > 1)
                writer.write("SampleId\t");

            writer.write("VariantCount\tGeneVariantCount\tRnaDepthSupport\tRnaDepthMedian\tRnaAlleleSupport\tRnaAlleleMedian");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private synchronized void writeSampleSummary(final String sampleId, final SampleVariantSummary summaryData)
    {
        try
        {
            if(mSampleIds.size() > 1)
                mWriter.write(format("%s\t", sampleId));

            int depthAboveThreshold = 0;
            int supportAboveThresholds = 0;
            for(VariantSupport support : summaryData.VariantFragments)
            {
                if(support.Depth >= mDepthThreshold)
                {
                    ++depthAboveThreshold;

                    if(support.AlleleFragments >= mSupportThreshold)
                        ++supportAboveThresholds;
                }
            }

            mWriter.write(format("%d\t%d\t%d\t%.1f\t%d\t%.1f",
                    summaryData.VariantCount, summaryData.GeneVariantCount,
                    depthAboveThreshold, summaryData.medianDepth(), supportAboveThresholds, summaryData.medianSupport()));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        addSampleIdFile(configBuilder, false);
        configBuilder.addPath(VCF_PATH, true, "Original VCF file");

        configBuilder.addInteger(MIN_DEPTH_FRAGS, "Required min RNA depth", DEFAULT_MIN_DEPTH_FRAGS);
        configBuilder.addInteger(MIN_SUPPORT_FRAGS, "Required min RNA allele fragments", DEFAULT_MIN_SUPPORT_FRAGS);

        configBuilder.addConfigItem(
                RNA_SAMPLE_SUFFIX, false,
                "RNA genotype sample name suffix, default: " + DEFAULT_RNA_SAMPLE_SUFFIX, DEFAULT_RNA_SAMPLE_SUFFIX);

        addThreadOptions(configBuilder);
        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        RnaHealthChecker rnaHealthChecker = new RnaHealthChecker(configBuilder);
        rnaHealthChecker.run();
    }
}

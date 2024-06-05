package com.hartwig.hmftools.esvee.utils;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.prep.PrepConfig.BLACKLIST_BED;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.esvee.prep.BlacklistLocations;

import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class BlacklistExplorer
{
    private final List<String> mSampleIds;
    private final String mVcfFilename;
    private final String mOutputFile;

    private final BlacklistLocations mBlacklistLocations;

    private final int mThreads;
    private final BufferedWriter mWriter;

    private static final String VCF_FILE = "vcf_file";
    private static final String OUTPUT_FILE = "output_file";

    public BlacklistExplorer(final ConfigBuilder configBuilder)
    {
        mSampleIds = loadSampleIdsFile(configBuilder);
        mVcfFilename = configBuilder.getValue(VCF_FILE);
        mBlacklistLocations = new BlacklistLocations(configBuilder.getValue(BLACKLIST_BED));
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);
        mWriter = initialiseWriter();
        mThreads = parseThreads(configBuilder);
    }

    public void run()
    {
        if(mVcfFilename == null)
        {
            SV_LOGGER.error("missing VCF");
            System.exit(1);
        }

        if(mSampleIds.isEmpty())
        {
            SV_LOGGER.error("no sample IDs specified");
            System.exit(1);
        }

        SV_LOGGER.info("searching for blacklist SVs in {} samples", mSampleIds.size());

        List<SampleTask> sampleTasks = Lists.newArrayList();
        for(int i = 0; i < min(mSampleIds.size(), mThreads); ++i)
        {
            sampleTasks.add(new SampleTask(i));
        }

        int taskIndex = 0;
        for(String sampleId : mSampleIds)
        {
            sampleTasks.get(taskIndex).sampleIds().add(sampleId);
            ++taskIndex;

            if(taskIndex >= sampleTasks.size())
                taskIndex = 0;
        }

        final List<Callable> callableList = sampleTasks.stream().collect(Collectors.toList());
        TaskExecutor.executeTasks(callableList, mThreads);

        closeBufferedWriter(mWriter);

        SV_LOGGER.info("blacklist SV search complete");
    }

    private synchronized void writeVariant(final String sampleId, final VariantContext variantContext, final BaseRegion blacklistRegion)
    {
        try
        {
            mWriter.write(format("%s,%s,%s,%d", sampleId, variantContext.getID(), variantContext.getContig(), variantContext.getStart()));

            String filtersStr = "";

            if(!variantContext.getFilters().isEmpty())
            {
                StringJoiner filters = new StringJoiner(";");
                variantContext.getFilters().forEach(x -> filters.add(x));
                filtersStr = filters.toString();
            }
            else
            {
                filtersStr = PASS;
            }

            int refDepth = variantContext.getAttributeAsInt(REF_DEPTH, 0)
                    + variantContext.getAttributeAsInt(REF_DEPTH_PAIR, 0);

            Genotype tumor = variantContext.getGenotype(1);

            int tumorFrags = getGenotypeAttributeAsInt(tumor, TOTAL_FRAGS, 0);

            mWriter.write(format(",%.1f,%d,%d,%s",
                    variantContext.getPhredScaledQual(), tumorFrags, refDepth, filtersStr));

            mWriter.write(format(",%d,%d", blacklistRegion.start(), blacklistRegion.end()));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            SV_LOGGER.info("writing comparison file: {}", mOutputFile);

            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            writer.write("SampleId,VcfId,Chromosome,Position");
            writer.write(",Qual,TumorFrags,RefDepth,Filters,RegionStart,RegionEnd");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
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

        public List<String> sampleIds() { return mSampleIds; }

        @Override
        public Long call()
        {
            SV_LOGGER.info("{}: processing {} samples", mTaskId, mSampleIds.size());

            int processed = 0;

            for(String sampleId : mSampleIds)
            {
                processSample(sampleId);

                ++processed;

                if((processed % 10) == 0)
                {
                    SV_LOGGER.debug("{}: processed {} samples", mTaskId, processed);
                }
            }

            SV_LOGGER.debug("{}: complete", mTaskId);

            return (long)0;
        }

        private void processSample(final String sampleId)
        {
            String vcfFilename = mVcfFilename.replaceAll("\\*", sampleId);

            if(!Files.exists(Paths.get(vcfFilename)))
            {
                SV_LOGGER.info("sample({}) vcf({}) doesn't exist, skipping", sampleId, vcfFilename);
                return;
            }

            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    vcfFilename, new VCFCodec(), false);

            try
            {
                int processed = 0;
                for(VariantContext variantContext : reader.iterator())
                {
                    String chromosome = variantContext.getContig();
                    int position = variantContext.getStart();

                    BaseRegion blacklistRegion = mBlacklistLocations.findBlacklistLocation(chromosome, position);

                    if(blacklistRegion != null)
                    {
                        writeVariant(sampleId, variantContext, blacklistRegion);
                    }

                    ++processed;

                    if((processed % 100000) == 0)
                    {
                        SV_LOGGER.debug("sample({}) processed {} variants", sampleId, processed);
                    }
                }
            }
            catch(IOException e)
            {
                SV_LOGGER.error("error reading vcf({}): {}", mVcfFilename, e.toString());
            }
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addSampleIdFile(configBuilder, false);
        configBuilder.addConfigItem(VCF_FILE, true, "VCF file, can use '*' in place of sampleIds");
        configBuilder.addPath(BLACKLIST_BED, false, "Blacklist BED file");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output filename");

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        BlacklistExplorer blacklistExplorer = new BlacklistExplorer(configBuilder);
        blacklistExplorer.run();
    }
}

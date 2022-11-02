package com.hartwig.hmftools.svprep.tools;

import static java.lang.Math.min;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READPAIR_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.REFERENCE_BREAKEND_READ_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKEND_COVERAGE;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.VARIANT_FRAGMENT_BREAKPOINT_COVERAGE;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConfig.BLACKLIST_BED;

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
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.svprep.BlacklistLocations;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
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

    public BlacklistExplorer(final CommandLine cmd)
    {
        mSampleIds = loadSampleIdsFile(cmd);
        mVcfFilename = cmd.getOptionValue(VCF_FILE);
        mBlacklistLocations = new BlacklistLocations(cmd.getOptionValue(BLACKLIST_BED));
        mOutputFile = cmd.getOptionValue(OUTPUT_FILE);
        mWriter = initialiseWriter();
        mThreads = parseThreads(cmd);
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

            int refDepth = variantContext.getAttributeAsInt(REFERENCE_BREAKEND_READ_COVERAGE, 0)
                    + variantContext.getAttributeAsInt(REFERENCE_BREAKEND_READPAIR_COVERAGE, 0);

            Genotype tumor = variantContext.getGenotype(1);
            int tumorFrags = getGenotypeAttributeAsInt(tumor, VARIANT_FRAGMENT_BREAKPOINT_COVERAGE, 0)
                    + getGenotypeAttributeAsInt(tumor, VARIANT_FRAGMENT_BREAKEND_COVERAGE, 0);

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

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addSampleIdFile(options);
        options.addOption(VCF_FILE, true, "VCF file, can use '*' in place of sampleIds");
        options.addOption(BLACKLIST_BED, true, "Blacklist BED file");
        options.addOption(OUTPUT_FILE, true, "Output filename");

        addOutputOptions(options);
        addLoggingOptions(options);
        addThreadOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        BlacklistExplorer blacklistExplorer = new BlacklistExplorer(cmd);
        blacklistExplorer.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}

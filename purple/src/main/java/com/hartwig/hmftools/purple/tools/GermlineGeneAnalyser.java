package com.hartwig.hmftools.purple.tools;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.utils.TaskExecutor;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class GermlineGeneAnalyser
{
    private final EnsemblDataCache mGeneDataCache;
    private final List<String> mSampleIds;
    private final String mPurpleDataDir;
    private final int mThreads;

    private final BufferedWriter mWriter;

    private static final String SAMPLE_ID_FILE = "sample_id_file";
    private static final String PURPLE_DATA_DIR = "purple_data_dir";
    private static final String THREADS = "threads";

    public GermlineGeneAnalyser(final CommandLine cmd)
    {
        mGeneDataCache = new EnsemblDataCache(cmd, V37);

        mGeneDataCache.setRequiredData(true, false, false, true);
        mGeneDataCache.load(false);

        mSampleIds = loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE));
        mPurpleDataDir = cmd.getOptionValue(PURPLE_DATA_DIR);
        mThreads = Integer.parseInt(cmd.getOptionValue(THREADS));

        mWriter = initialiseWriter(cmd.getOptionValue(OUTPUT_DIR));
    }

    public void run()
    {
        if(mSampleIds.isEmpty())
        {
            PPL_LOGGER.error("missing sampleIds, exiting");
            System.exit(1);
        }

        PPL_LOGGER.info("running Purple germline gene analysis for {} samples", mSampleIds.size());

        List<SampleGermlineGeneTask> sampleTasks = Lists.newArrayList();

        if(mThreads > 1)
        {
            for(int i = 0; i < min(mSampleIds.size(), mThreads); ++i)
            {
                sampleTasks.add(new SampleGermlineGeneTask(i, mWriter, mGeneDataCache, mPurpleDataDir));
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
            SampleGermlineGeneTask sampleTask = new SampleGermlineGeneTask(0, mWriter, mGeneDataCache, mPurpleDataDir);

            sampleTask.getSampleIds().addAll(mSampleIds);

            sampleTasks.add(sampleTask);

            sampleTask.call();
        }

        closeBufferedWriter(mWriter);

        PPL_LOGGER.info("Purple germline gene analysis cnmplete");
    }

    private BufferedWriter initialiseWriter(final String outputDir)
    {
        try
        {
            String fileName = outputDir + "purple_germline_gene_data.csv";

            // PPL_LOGGER.info("writing output file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId,Chromosome,RegionStart,RegionEnd,GermlineStatus,RegionBafCount,RegionDWC");
            writer.write(",RegionObsTumorRatio,RegionObsNormalRatio,RegionRefNormalisedCN");
            writer.write(",PurpleCN,PurpleDWC,GcContent,MajorAlleleCN,MinorAlleleCN");
            writer.write(",GeneName,TransName,TransStart,TransEnd,ExonCount,ExonRankStart,ExonRankEnd");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to initialise output file output: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeGeneOverlapData(
            final BufferedWriter writer, final String sampleId, final ObservedRegion region, final PurpleCopyNumber copyNumber,
            final GeneData geneData, final TranscriptData transData, final List<ExonData> overlappedExons)
    {
        try
        {
            writer.write(String.format("%s,%s,%d,%d,%s,%d,%d,%.3f,%.3f,%.3f",
                    sampleId, region.chromosome(), region.start(), region.end(), region.germlineStatus(),
                    region.bafCount(), region.depthWindowCount(), region.observedTumorRatio(), region.observedNormalRatio(),
                    region.refNormalisedCopyNumber()));

            writer.write(String.format(",%.2f,%d,%.2f,%.2f,%.2f",
                    copyNumber.averageTumorCopyNumber(), copyNumber.depthWindowCount(), copyNumber.gcContent(),
                    copyNumber.majorAlleleCopyNumber(), copyNumber.minorAlleleCopyNumber()));

            int exonRangeMin = overlappedExons.stream().mapToInt(x -> x.Rank).min().orElse(0);
            int exonRangeMax = overlappedExons.stream().mapToInt(x -> x.Rank).max().orElse(0);

            writer.write(String.format(",%s,%s,%d,%d,%d,%d,%d",
                    geneData.GeneName, transData.TransName, transData.TransStart, transData.TransEnd,
                    transData.exons().size(), exonRangeMin, exonRangeMax));

            writer.newLine();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write germline gene overlap file output: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE_ID_FILE, true, "Sample ID file");
        options.addOption(THREADS, true, "Thread count, default = 0 (disabled)");

        options.addOption(PURPLE_DATA_DIR, true, "Directory pattern for sample purple directory");
        addEnsemblDir(options);
        addLoggingOptions(options);
        addOutputOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        GermlineGeneAnalyser germlineGeneAnalyser = new GermlineGeneAnalyser(cmd);
        germlineGeneAnalyser.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
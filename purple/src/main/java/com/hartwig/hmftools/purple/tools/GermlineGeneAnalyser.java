package com.hartwig.hmftools.purple.tools;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeVersion;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addSampleIdFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
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
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.common.perf.TaskExecutor;

import org.jetbrains.annotations.NotNull;

public class GermlineGeneAnalyser
{
    private final EnsemblDataCache mGeneDataCache;
    private final List<String> mSampleIds;
    private final String mPurpleDataDir;
    private final int mThreads;

    private final BufferedWriter mWriter;
    private final boolean mWriteVerbose;

    private static final String WRITE_VERBOSE = "write_verbose";

    public GermlineGeneAnalyser(final ConfigBuilder configBuilder)
    {
        mGeneDataCache = new EnsemblDataCache(configBuilder);

        mGeneDataCache.setRequiredData(true, false, false, true);
        mGeneDataCache.load(false);

        mSampleIds = loadSampleIdsFile(configBuilder);
        mPurpleDataDir = configBuilder.getValue(PURPLE_DIR_CFG);
        mThreads = parseThreads(configBuilder);

        mWriteVerbose = configBuilder.hasFlag(WRITE_VERBOSE);
        mWriter = initialiseWriter(parseOutputDir(configBuilder));
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
                sampleTasks.add(new SampleGermlineGeneTask(i, mWriter, mGeneDataCache, mPurpleDataDir, mWriteVerbose));
            }

            int taskIndex = 0;
            for(String sampleId : mSampleIds)
            {
                if(taskIndex >= sampleTasks.size())
                    taskIndex = 0;

                sampleTasks.get(taskIndex).getSampleIds().add(sampleId);

                ++taskIndex;
            }

            final List<Callable<Void>> callableList = sampleTasks.stream().collect(Collectors.toList());
            TaskExecutor.executeTasks(callableList, mThreads);
        }
        else
        {
            SampleGermlineGeneTask sampleTask = new SampleGermlineGeneTask(0, mWriter, mGeneDataCache, mPurpleDataDir, mWriteVerbose);

            sampleTask.getSampleIds().addAll(mSampleIds);

            sampleTasks.add(sampleTask);

            sampleTask.call();
        }

        closeBufferedWriter(mWriter);

        PPL_LOGGER.info("Purple germline gene analysis complete");
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

            if(mWriteVerbose)
            {
                writer.write(",PurpleCN,PurpleDWC,GcContent,MajorAlleleCN,MinorAlleleCN");
                writer.write(",GeneName,TransName,TransStart,TransEnd,ExonCount,ExonRankStart,ExonRankEnd");
            }
            else
            {
                writer.write(",GeneNames");
            }

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to initialise output file output: {}", e.toString());
            return null;
        }
    }

    public synchronized static void writeGeneDeletionData(
            final BufferedWriter writer, final String sampleId, final ObservedRegion region, final String geneNames)
    {
        try
        {
            writer.write(String.format("%s,%s,%d,%d,%s,%d,%d,%.3f,%.3f,%.3f,%s",
                    sampleId, region.chromosome(), region.start(), region.end(), region.germlineStatus(),
                    region.bafCount(), region.depthWindowCount(), region.observedTumorRatio(), region.observedNormalRatio(),
                    region.refNormalisedCopyNumber(), geneNames));

            writer.newLine();
        }
        catch(IOException e)
        {
            PPL_LOGGER.error("failed to write germline gene deletion file output: {}", e.toString());
        }
    }

    public synchronized static void writeGeneDeletionDetails(
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

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        addSampleIdFile(configBuilder, true);

        configBuilder.addPath(PURPLE_DIR_CFG, true, "Directory pattern for sample purple directory");
        configBuilder.addFlag(WRITE_VERBOSE, "Write transcript and exon details for each deleted gene segment");
        addRefGenomeVersion(configBuilder);

        addEnsemblDir(configBuilder);
        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
        addThreadOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        setLogLevel(configBuilder);

        GermlineGeneAnalyser germlineGeneAnalyser = new GermlineGeneAnalyser(configBuilder);
        germlineGeneAnalyser.run();
    }
}
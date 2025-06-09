package com.hartwig.hmftools.geneutils.utils;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.region.PartitionUtils.buildPartitions;
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;
import org.jetbrains.annotations.NotNull;

public class GenomeMappabilityFinder
{
    private final String mOutputFile;

    private final BufferedWriter mWriter;
    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenVersion;

    private final int mBaseWindowLength;
    private final int mBaseBufferLength;

    private final int mPartitionSize;
    private final int mThreads;
    private final SpecificRegions mSpecificChrRegions;

    private final BwaMemAligner mAligner;

    // config items
    private static final String CFG_OUTPUT_FILE = "output_file";
    protected static final String CFG_BASE_WINDOW_LENGTH = "window_length";
    protected static final String CFG_BASE_BUFFER_LENGTH = "buffer_length";
    protected static final String CFG_PARTITION_SIZE = "partition_size";

    private static final int DEFAULT_BASE_LENGTH = 20;
    private static final int DEFAULT_BUFFER_LENGTH = 40;

    public GenomeMappabilityFinder(final ConfigBuilder configBuilder)
    {
        mOutputFile = configBuilder.getValue(CFG_OUTPUT_FILE);
        mWriter = initialiseWriter();

        mRefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));
        mRefGenVersion = deriveRefGenomeVersion(mRefGenome);

        mPartitionSize = configBuilder.getInteger(CFG_PARTITION_SIZE);
        mBaseWindowLength = configBuilder.getInteger(CFG_BASE_WINDOW_LENGTH);
        mBaseBufferLength = configBuilder.getInteger(CFG_BASE_BUFFER_LENGTH);

        mThreads = parseThreads(configBuilder);
        mSpecificChrRegions = SpecificRegions.from(configBuilder);

        String refGenomeImageFile = configBuilder.getValue(REF_GENOME) + ".img";

        loadAlignerLibrary(configBuilder.getValue(BWA_LIB_PATH));
        mAligner = initialiseBwaAligner(refGenomeImageFile);
    }

    private BwaMemAligner initialiseBwaAligner(final String refGenomeImageFile)
    {
        if(refGenomeImageFile.isEmpty() || !Files.exists(Paths.get(refGenomeImageFile)))
            System.exit(1);

        try
        {
            BwaMemIndex index = new BwaMemIndex(refGenomeImageFile);
            BwaMemAligner aligner = new BwaMemAligner(index);

            // adjust parameters as required

            GU_LOGGER.debug("BWA align options: ");
            GU_LOGGER.debug("  BWA Bandwidth: {}", aligner.getBandwidthOption());
            GU_LOGGER.debug("  BWA MatchScore: {}", aligner.getMatchScoreOption());
            GU_LOGGER.debug("  BWA OutputScoreThreshold: {}", aligner.getOutputScoreThresholdOption());
            GU_LOGGER.debug("  BWA SplitFactor: {}", aligner.getSplitFactorOption());
            GU_LOGGER.debug("  BWA Mismatch: {}", aligner.getXADropRatio());
            GU_LOGGER.debug("  BWA IGapOpen: {}", aligner.getIGapOpenPenaltyOption());
            GU_LOGGER.debug("  BWA IGapOpenExtend: {}", aligner.getIGapExtendPenaltyOption());
            GU_LOGGER.debug("  BWA DGapOpen: {}", aligner.getDGapOpenPenaltyOption());
            GU_LOGGER.debug("  BWA DGapOpenExtend: {}", aligner.getDGapExtendPenaltyOption());
            GU_LOGGER.debug("  BWA DropRatio: {}", aligner.getDropRatioOption());
            GU_LOGGER.debug("  BWA ZDrop: {}", aligner.getZDropOption());
            GU_LOGGER.debug("  BWA XADropRatio: {}", aligner.getXADropRatio());
            GU_LOGGER.debug("  BWA MaxXAHits: {}", aligner.getMaxXAHitsOption());
            GU_LOGGER.debug("  BWA MaxXAHitsAlt: {}", aligner.getMaxXAHitsAltOption());

            // aligner.setMismatchPenaltyOption(5);

            return aligner;
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to initialise BWA aligner: {}", e.toString());
            return null;
        }
    }
    public void run()
    {
        GU_LOGGER.info("analysing blacklist regions");

        long startTimeMs = System.currentTimeMillis();

        RefGenomeCoordinates coordinates = mRefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        List<ChrBaseRegion> regions = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            String chrStr = mRefGenVersion.versionedChromosome(chromosome.toString());

            if(mSpecificChrRegions.excludeChromosome(chrStr))
                continue;

            List<ChrBaseRegion> chrRegions = buildPartitions(chrStr, coordinates.length(chrStr), mPartitionSize);

            for(ChrBaseRegion region : chrRegions)
            {
                if(mSpecificChrRegions.hasFilters() && !mSpecificChrRegions.includeRegion(region))
                    continue;

                regions.add(region);
            }
        }

        List<RegionTask> regionTasks = regions.stream().map(x -> new RegionTask(x)).collect(Collectors.toList());

        List<Callable> callableList = regionTasks.stream().collect(Collectors.toList());
        if(!TaskExecutor.executeTasks(callableList, mThreads))
            System.exit(1);

        closeBufferedWriter(mWriter);

        GU_LOGGER.info("Genome mappability analysis complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private class RegionTask implements Callable
    {
        private final ChrBaseRegion mRegion;

        public RegionTask(final ChrBaseRegion region)
        {
            mRegion = region;
        }

        @Override
        public Long call()
        {
            int currentWindowStart = mRegion.start();
            int windowEnd = currentWindowStart + mBaseWindowLength - 1;

            int minRegionEnd = mRegion.end() - (int)(mBaseWindowLength * 0.1);

            while(windowEnd <= minRegionEnd)
            {
                analyseWindow(new BaseRegion(currentWindowStart, min(windowEnd, mRegion.end())));

                currentWindowStart = windowEnd + 1;
                windowEnd = currentWindowStart + mBaseWindowLength - 1;
            }

            return (long)0;
        }

        private void analyseWindow(final BaseRegion region)
        {
            byte[] refBases = mRefGenome.getBases(mRegion.Chromosome, region.start(), region.end());

            // writeWindowData(mRegion.Chromosome);
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(FLD_CHROMOSOME).add(FLD_POSITION_START).add(FLD_POSITION_END);

            // sj.add("ScoreInfo") etc;

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private synchronized void writeWindowData(
            final String chromosome, final int positionStart, final int positionEnd)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(chromosome).add(String.valueOf(positionStart)).add(String.valueOf(positionEnd));

            mWriter.write(sj.toString());
            mWriter.newLine();

        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addRefGenomeFile(configBuilder, true);

        configBuilder.addConfigItem(CFG_OUTPUT_FILE, true, "Output filename");
        configBuilder.addInteger(CFG_PARTITION_SIZE, "Read length", 1_000_000);

        configBuilder.addInteger(CFG_BASE_WINDOW_LENGTH, "Base window length for analysis", DEFAULT_BASE_LENGTH);
        configBuilder.addInteger(CFG_BASE_BUFFER_LENGTH, "Base window length for analysis", DEFAULT_BUFFER_LENGTH);

        addOutputOptions(configBuilder);
        addLoggingOptions(configBuilder);
        addThreadOptions(configBuilder);
        addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        GenomeMappabilityFinder genomeMappabilityFinder = new GenomeMappabilityFinder(configBuilder);
        genomeMappabilityFinder.run();
    }
}

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
import static com.hartwig.hmftools.common.region.SpecificRegions.addSpecificChromosomesRegionsConfig;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

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

// TODO: doc
public class OffTargetRiskProfiler
{
    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenVersion;
    // Genome region filtering.
    private final SpecificRegions mSpecificChrRegions;

    // Genome is partitioned in windows of this length, each of which are given a score.
    private final int mBaseWindowLength;
    // Each window is spaced this far apart. < mBaseWindowLength implies overlapping windows.
    private final int mBaseWindowSpacing;

    // Alignment score must exceed this to count towards the risk score.
    private final int mMatchScoreThreshold;
    // Amount that 1 alignment match counts towards the risk score.
    // E.g. value of 10 means an alignment with score = mMatchScoreThreshold contributes 10 risk score points.
    private final int mMatchScoreOffset;

    // Performance tuning parameters.
    // How many windows to align with BWA-MEM at once.
    private final int mBatchSize;
    private final int mThreads;

    private final BwaMemAligner mAligner;

    private final BufferedWriter mOutputWriter;

    private static final String CFG_OUTPUT_FILE = "output_file";

    private static final String CFG_BASE_WINDOW_LENGTH = "window_length";
    private static final int DEFAULT_WINDOW_LENGTH = 40;

    private static final String CFG_BASE_WINDOW_SPACING = "window_spacing";
    private static final int DEFAULT_WINDOW_SPACING = 40;

    private static final String CFG_MATCH_SCORE_THRESHOLD = "match_score_threshold";
    private static final int DEFAULT_MATCH_SCORE_THRESHOLD = 24;    // Determined empirically via experiment

    private static final String CFG_MATCH_SCORE_OFFSET = "match_score_offset";
    private static final int DEFAULT_MATCH_SCORE_OFFSET = 26;   // Determined empirically via experiment

    private static final String CFG_BATCH_SIZE = "batch_size";
    private static final int DEFAULT_BATCH_SIZE = 100000;

    private static final String FLD_RISK_SCORE = "risk_score";
    private static final String FLD_QUALITY_SCORE = "quality_score";

    public OffTargetRiskProfiler(final ConfigBuilder configBuilder)
    {
        mRefGenome = loadRefGenome(configBuilder.getValue(REF_GENOME));
        mRefGenVersion = deriveRefGenomeVersion(mRefGenome);

        mSpecificChrRegions = SpecificRegions.from(configBuilder);

        mBaseWindowLength = configBuilder.getInteger(CFG_BASE_WINDOW_LENGTH);
        mBaseWindowSpacing = configBuilder.getInteger(CFG_BASE_WINDOW_SPACING);
        mMatchScoreThreshold = configBuilder.getInteger(CFG_MATCH_SCORE_THRESHOLD);
        mMatchScoreOffset = configBuilder.getInteger(CFG_MATCH_SCORE_OFFSET);

        mBatchSize = configBuilder.getInteger(CFG_BATCH_SIZE);
        mThreads = parseThreads(configBuilder);

        String refGenomeImageFile = configBuilder.getValue(REF_GENOME) + ".img";

        loadAlignerLibrary(configBuilder.getValue(BWA_LIB_PATH));
        mAligner = initialiseBwaAligner(refGenomeImageFile);

        mOutputWriter = initialiseOutputWriter(configBuilder.getValue(CFG_OUTPUT_FILE));
    }

    private BwaMemAligner initialiseBwaAligner(final String refGenomeImageFile)
    {
        if(refGenomeImageFile.isEmpty() || !Files.exists(Paths.get(refGenomeImageFile)))
            System.exit(1);

        try
        {
            BwaMemIndex index = new BwaMemIndex(refGenomeImageFile);
            BwaMemAligner aligner = new BwaMemAligner(index);

            // Output many alignments per query
            aligner.setFlagOption(aligner.getFlagOption() | aligner.MEM_F_ALL);
            aligner.setOutputScoreThresholdOption(20);
            // Don't prune seeds with many occurrences in the genome. This is a key performance tuning parameter.
            aligner.setMaxMemIntvOption(2000);
            aligner.setMaxSeedOccurencesOption(10000);
            // Other minor params to encourage more alignments to be found.
            aligner.setBandwidthOption(mBaseWindowLength);
            aligner.setZDropOption(500);
            aligner.setSplitFactorOption(0.5f);
            aligner.setDropRatioOption(0.1f);

            GU_LOGGER.debug("BWA-MEM options:");
            GU_LOGGER.debug("  MinSeedLength: {}", aligner.getMinSeedLengthOption());
            GU_LOGGER.debug("  SplitFactor: {}", aligner.getSplitFactorOption());
            GU_LOGGER.debug("  SplitWidth: {}", aligner.getSplitWidthOption());
            GU_LOGGER.debug("  MaxSeedOccurrences: {}", aligner.getMaxSeedOccurencesOption());
            GU_LOGGER.debug("  MaxMemOccurrences: {}", aligner.getMaxMemIntvOption());
            GU_LOGGER.debug("  DropRatio: {}", aligner.getDropRatioOption());
            GU_LOGGER.debug("  Match: {}", aligner.getMatchScoreOption());
            GU_LOGGER.debug("  Mismatch: {}", aligner.getXADropRatio());
            GU_LOGGER.debug("  IGapOpen: {}", aligner.getIGapOpenPenaltyOption());
            GU_LOGGER.debug("  IGapExtend: {}", aligner.getIGapExtendPenaltyOption());
            GU_LOGGER.debug("  DGapOpen: {}", aligner.getDGapOpenPenaltyOption());
            GU_LOGGER.debug("  DGapExtend: {}", aligner.getDGapExtendPenaltyOption());
            GU_LOGGER.debug("  Clip3: {}", aligner.getClip3PenaltyOption());
            GU_LOGGER.debug("  Clip5: {}", aligner.getClip5PenaltyOption());
            GU_LOGGER.debug("  Bandwidth: {}", aligner.getBandwidthOption());
            GU_LOGGER.debug("  ZDrop: {}", aligner.getZDropOption());
            GU_LOGGER.debug("  OutputScoreThreshold: {}", aligner.getOutputScoreThresholdOption());
            GU_LOGGER.debug("  Flags: {}", aligner.getFlagOption());

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

        Stream<ChrBaseRegion> regions = getWindowRegions();

        // TODO
        List<RegionTask> regionTasks = regions.stream().map(x -> new RegionTask(x)).collect(Collectors.toList());
        List<Callable> callableList = regionTasks.stream().collect(Collectors.toList());
        if(!TaskExecutor.executeTasks(callableList, mThreads))
            System.exit(1);

        closeBufferedWriter(mOutputWriter);

        GU_LOGGER.info("Genome mappability analysis complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private Stream<ChrBaseRegion> getWindowRegions() {
        RefGenomeCoordinates coordinates = mRefGenVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;

        return Arrays.stream(HumanChromosome.values())
                .map(chr -> mRefGenVersion.versionedChromosome(chr.toString()))
                .filter(mSpecificChrRegions::includeChromosome)
                .flatMap(chr -> createWindowRegions(chr, coordinates.length(chr)))
                .filter(region -> !mSpecificChrRegions.hasFilters() || mSpecificChrRegions.includeRegion(region));
    }

    // Partition a chromosome into windows to analyse.
    private Stream<ChrBaseRegion> createWindowRegions(String chromosome, int chromosomeLength) {
        // TODO: is start and end 0 or 1 indexed? inclusive/exclusive?
        return IntStream.iterate(1, start -> true, start -> start + mBaseWindowSpacing)
                .mapToObj(start -> new ChrBaseRegion(chromosome, start, start + mBaseWindowLength - 1))
                .takeWhile(region -> region.end() <= chromosomeLength);
    }

    // TODO
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

    private BufferedWriter initialiseOutputWriter(String path)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(path, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME).add(FLD_POSITION_START).add(FLD_POSITION_END)
                    .add(FLD_RISK_SCORE).add(FLD_QUALITY_SCORE);
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

    private void writeWindowData(
            final String chromosome, final int positionStart, final int positionEnd,
            final double riskScore, final double qualityScore)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(chromosome).add(String.valueOf(positionStart)).add(String.valueOf(positionEnd))
                    .add(String.valueOf(riskScore)).add(String.valueOf(qualityScore));
            mOutputWriter.write(sj.toString());
            mOutputWriter.newLine();
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
        addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.addInteger(CFG_BASE_WINDOW_LENGTH, "Base window length for analysis", DEFAULT_WINDOW_LENGTH);
        configBuilder.addInteger(CFG_BASE_WINDOW_SPACING, "Offset through the genome of each base window", DEFAULT_WINDOW_SPACING);
        configBuilder.addInteger(CFG_MATCH_SCORE_THRESHOLD, "Minimum alignment score to consider a match against the window", DEFAULT_MATCH_SCORE_THRESHOLD);
        configBuilder.addInteger(CFG_MATCH_SCORE_OFFSET, "Points contributed when alignment score = threshold", DEFAULT_MATCH_SCORE_OFFSET);

        configBuilder.addInteger(CFG_BATCH_SIZE, "Number of windows to align simultaneously", DEFAULT_BATCH_SIZE);
        addThreadOptions(configBuilder);

        configBuilder.addConfigItem(CFG_OUTPUT_FILE, true, "Output filename");

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        OffTargetRiskProfiler offTargetRiskProfiler = new OffTargetRiskProfiler(configBuilder);
        offTargetRiskProfiler.run();
    }
}

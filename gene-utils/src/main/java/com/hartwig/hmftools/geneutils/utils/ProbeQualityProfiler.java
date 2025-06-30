package com.hartwig.hmftools.geneutils.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH_DESC;
import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
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
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

// Tool for helping with panel probe design. It produces a file which annotates the genome with information informing how likely a probe is
// to hybridise with off-target genome regions.
// The genome is partitioned into small base windows which are scored according to a probe off-target risk model.
// The model uses BWA-MEM alignments as a heuristic for the chance of off-target hybridisation, based on the assumption that hybridisation
// occurs with a sufficiently long exact match. The model was developed and validated against data from a previous panel design.
// Output is a TSV file containing scores for each window in the genome.
public class ProbeQualityProfiler
{
    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenomeVersion;
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
    // If true, output more information than just the final quality score. Useful for debugging.
    private final boolean mVerboseOutput;
    // Number of decimals included in the quality score.
    private static final int QUALITY_SCORE_PRECISION = 4;

    private static final String CFG_BASE_WINDOW_LENGTH = "window_length";
    private static final int DEFAULT_WINDOW_LENGTH = 40;

    private static final String CFG_BASE_WINDOW_SPACING = "window_spacing";
    private static final int DEFAULT_WINDOW_SPACING = 20;

    private static final String CFG_MATCH_SCORE_THRESHOLD = "match_score_threshold";
    private static final int DEFAULT_MATCH_SCORE_THRESHOLD = 24;    // Determined empirically via experiment

    private static final String CFG_MATCH_SCORE_OFFSET = "match_score_offset";
    private static final int DEFAULT_MATCH_SCORE_OFFSET = 26;   // Determined empirically via experiment

    private static final String CFG_BATCH_SIZE = "batch_size";
    private static final int DEFAULT_BATCH_SIZE = 100000;

    private static final String CFG_OUTPUT_FILE = "output_file";

    private static final String CFG_VERBOSE_OUTPUT = "verbose_output";

    private static final String FLD_QUALITY_SCORE = "QualityScore";
    private static final String FLD_RISK_SCORE = "RiskScore";
    private static final String FLD_OFF_TARGET_COUNT = "OffTargetCount";
    private static final String FLD_OFF_TARGET_SCORE_SUM = "OffTargetScoreSum";

    public ProbeQualityProfiler(final ConfigBuilder configBuilder)
    {
        String refGenome = configBuilder.getValue(REF_GENOME);
        GU_LOGGER.debug("Ref genome: {}", refGenome);
        mRefGenome = loadRefGenome(refGenome);
        mRefGenomeVersion = deriveRefGenomeVersion(mRefGenome);
        GU_LOGGER.debug("Ref genome version: {}", mRefGenomeVersion);

        mSpecificChrRegions = SpecificRegions.from(configBuilder);
        if (mSpecificChrRegions == null) {
            throw new RuntimeException("Invalid config");
        }

        mBaseWindowLength = configBuilder.getInteger(CFG_BASE_WINDOW_LENGTH);
        if (mBaseWindowLength < 15) {
            // Less than 15 bases probably doesn't make much sense and will cause issues trying to find appropriate params for BWA-MEM.
            throw new RuntimeException(String.format("%s must be >= 10", CFG_BASE_WINDOW_LENGTH));
        }
        GU_LOGGER.debug("Base window length: {}", mBaseWindowLength);

        mBaseWindowSpacing = configBuilder.getInteger(CFG_BASE_WINDOW_SPACING);
        if (mBaseWindowSpacing < 1) {
            throw new RuntimeException(String.format("%s must be >= 1", CFG_BASE_WINDOW_SPACING));
        }
        GU_LOGGER.debug("Base window spacing: {}", mBaseWindowSpacing);

        mMatchScoreThreshold = configBuilder.getInteger(CFG_MATCH_SCORE_THRESHOLD);
        if (mMatchScoreThreshold > mBaseWindowLength) {
            // If this is true then all alignments will be excluded which is useless.
            throw new RuntimeException(String.format("%s must be <= %s", CFG_MATCH_SCORE_THRESHOLD, CFG_BASE_WINDOW_LENGTH));
        }
        GU_LOGGER.debug("Match score threshold: {}", mMatchScoreThreshold);

        mMatchScoreOffset = configBuilder.getInteger(CFG_MATCH_SCORE_OFFSET);
        if (mMatchScoreOffset < 0) {
            // Negative values will break the risk model maths.
            throw new RuntimeException(String.format("%s must be >= 0", CFG_MATCH_SCORE_OFFSET));
        }
        GU_LOGGER.debug("Match score offset: {}", mMatchScoreOffset);

        mBatchSize = configBuilder.getInteger(CFG_BATCH_SIZE);
        if (mBatchSize < 1) {
            throw new RuntimeException(String.format("%s must be >= 1", CFG_BATCH_SIZE));
        }
        GU_LOGGER.debug("Batch size: {}", mBatchSize);

        mThreads = TaskExecutor.parseThreads(configBuilder);
        if (mThreads < 1) {
            throw new RuntimeException(String.format("%s must be >= 1", TaskExecutor.THREADS));
        }
        GU_LOGGER.debug("Threads: {}", mThreads);

        mVerboseOutput = configBuilder.hasFlag(CFG_VERBOSE_OUTPUT);
        GU_LOGGER.debug("Verbose output: {}", mVerboseOutput);

        String refGenomeImageFile = refGenome + ".img";

        loadAlignerLibrary(configBuilder.getValue(BWA_LIB_PATH));
        mAligner = initialiseBwaAligner(refGenomeImageFile);

        String outputFile = configBuilder.getValue(CFG_OUTPUT_FILE);
        GU_LOGGER.debug("Output file: {}", outputFile);
        mOutputWriter = initialiseOutputWriter(outputFile);
    }

    private BwaMemAligner initialiseBwaAligner(final String refGenomeImageFile)
    {
        if(refGenomeImageFile.isEmpty() || !Files.exists(Paths.get(refGenomeImageFile)))
            System.exit(1);

        try
        {
            BwaMemIndex index = new BwaMemIndex(refGenomeImageFile);
            BwaMemAligner aligner = new BwaMemAligner(index);

            // Ensure we can find alignments fitting our parameters.
            aligner.setMinSeedLengthOption(min(min(19, mMatchScoreThreshold + 10), mBaseWindowLength / 2));
            // Output many alignments per query.
            aligner.setFlagOption(aligner.getFlagOption() | BwaMemAligner.MEM_F_ALL);
            aligner.setOutputScoreThresholdOption(mMatchScoreThreshold);
            // Don't prune seeds with many occurrences in the genome. This is a key performance tuning parameter.
            aligner.setMaxMemIntvOption(2000);
            aligner.setMaxSeedOccurencesOption(10000);
            // Other minor params to encourage more alignments to be found.
            aligner.setBandwidthOption(mBaseWindowLength);
            aligner.setSplitFactorOption(0.5f);
            aligner.setDropRatioOption(0.1f);
            // Performance params.
            aligner.setNThreadsOption(mThreads);

            GU_LOGGER.debug("BWA-MEM options:");
            GU_LOGGER.debug("  MinSeedLength: {}", aligner.getMinSeedLengthOption());
            GU_LOGGER.debug("  SplitFactor: {}", aligner.getSplitFactorOption());
            GU_LOGGER.debug("  SplitWidth: {}", aligner.getSplitWidthOption());
            GU_LOGGER.debug("  MaxSeedOccurrences: {}", aligner.getMaxSeedOccurencesOption());
            GU_LOGGER.debug("  MaxMemOccurrences: {}", aligner.getMaxMemIntvOption());
            GU_LOGGER.debug("  DropRatio: {}", aligner.getDropRatioOption());
            GU_LOGGER.debug("  Match: {}", aligner.getMatchScoreOption());
            GU_LOGGER.debug("  Mismatch: {}", aligner.getMismatchPenaltyOption());
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
            GU_LOGGER.debug("  Threads: {}", aligner.getNThreadsOption());

            return aligner;
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to initialise BWA aligner: {}", e.toString());
            return null;
        }
    }

    private BufferedWriter initialiseOutputWriter(String path)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(path, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME).add(FLD_POSITION_START).add(FLD_POSITION_END)
                    .add(FLD_QUALITY_SCORE);
            if (mVerboseOutput) {
                sj.add(FLD_RISK_SCORE);
                sj.add(FLD_OFF_TARGET_COUNT);
                sj.add(FLD_OFF_TARGET_SCORE_SUM);
            }
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

    public void run()
    {
        GU_LOGGER.info("Starting");

        long startTimeMs = System.currentTimeMillis();

        Stream<ChrBaseRegion> regions = createBaseWindowRegions();
        ProcessingStats stats = processBaseWindows(regions);

        closeBufferedWriter(mOutputWriter);

        GU_LOGGER.info("Analysis complete, mins({})", runTimeMinsStr(startTimeMs));
        GU_LOGGER.info("Window stats:");
        GU_LOGGER.info("  Total: {}", stats.totalWindows);
        GU_LOGGER.info("  Analysed: {}", stats.totalWindows);
        GU_LOGGER.info("  Denormal: {}", stats.denormalWindows);
    }

    private Stream<ChrBaseRegion> createBaseWindowRegions() {
        GU_LOGGER.info("Creating base window region stream");
        RefGenomeCoordinates coordinates = mRefGenomeVersion.is37() ? RefGenomeCoordinates.COORDS_37 : RefGenomeCoordinates.COORDS_38;
        return Arrays.stream(HumanChromosome.values())
                .map(chr -> mRefGenomeVersion.versionedChromosome(chr.toString()))
                .filter(mSpecificChrRegions::includeChromosome)
                .flatMap(chr -> createBaseWindowRegions(chr, coordinates.length(chr)));
    }

    private Stream<ChrBaseRegion> createBaseWindowRegions(String chromosome, int chromosomeLength) {
        GU_LOGGER.debug("Creating base window stream for chromosome: {}", chromosome);
        // Generate the windows from the configured specific regions, rather than enumerating all windows and filtering, for performance.
        List<ChrBaseRegion> regions = mSpecificChrRegions.Regions.stream()
                .filter(region -> region.chromosome().equals(chromosome)).collect(Collectors.toList());
        if (regions.isEmpty()) {
            regions.add(new ChrBaseRegion(chromosome, 1, chromosomeLength - 1));
        }
        return regions.stream()
                .flatMap(this::createBaseWindowRegions)
                // Might be a better way to handle the end of chromosomes, for now exclude them if the window doesn't line up.
                .takeWhile(region -> region.end() <= chromosomeLength);
    }

    // Create base window regions which fully cover the specified region.
    private Stream<ChrBaseRegion> createBaseWindowRegions(ChrBaseRegion region) {
        GU_LOGGER.debug("Creating base window stream for region: {}", region);
        // First window may start before the start of the specified region.
        int initial = baseWindowStartCoveringPosition(region.start());
        // Last window could extend past the end of the specified region.
        return IntStream.iterate(initial, start -> start <= region.end(), start -> start + mBaseWindowSpacing)
                .mapToObj(start -> new ChrBaseRegion(region.chromosome(), start, start + mBaseWindowLength - 1));
    }

    // Finds the start position of the nearest base window which covers the given position.
    private int baseWindowStartCoveringPosition(int position) {
        int position0Idx = position - 1;
        int mod = position0Idx % mBaseWindowSpacing;
        return max(position - mod, 1);
    }

    private record ProcessingStats(
            // Total number of windows considered.
            long totalWindows,
            // Number of windows analysed (i.e. not skipped).
            long analysedWindows,
            // How many base windows were skipped because they contained bases which are not ACGT.
            long denormalWindows
    ) {
        ProcessingStats() {
            this(0, 0, 0);
        }

        ProcessingStats add(ProcessingStats other) {
            return new ProcessingStats(
                    totalWindows + other.totalWindows,
                    analysedWindows + other.analysedWindows,
                    denormalWindows + other.denormalWindows
            );
        }
    }

    private ProcessingStats processBaseWindows(Stream<ChrBaseRegion> regions) {
        GU_LOGGER.info("Processing base windows");
        return partitionRegionsIntoBatches(regions)
                .map(this::processBaseWindowBatch)
                .reduce(new ProcessingStats(), ProcessingStats::add);
    }

    private Stream<List<ChrBaseRegion>> partitionRegionsIntoBatches(Stream<ChrBaseRegion> regions) {
        Iterator<ChrBaseRegion> iterator = regions.iterator();
        return Stream.generate(() -> {
            List<ChrBaseRegion> batch = new ArrayList<>(mBatchSize);
            for (int i = 0; i < mBatchSize && iterator.hasNext(); i++) {
                batch.add(iterator.next());
            }
            return batch;
        }).takeWhile(b -> !b.isEmpty());
    }

    private ProcessingStats processBaseWindowBatch(List<ChrBaseRegion> batch) {
        GU_LOGGER.debug("Processing base window batch of size {}", batch.size());
        GU_LOGGER.debug("First base window: {}", batch.get(0).toString());

        GU_LOGGER.debug("Retrieving base window sequences");
        // Somewhat awkward stream processing here because we need the base sequence to filter on
        // and also need to keep the region info for later.
        List<RegionWithSequence> preprocessedRegions = batch.stream()
                .map(region -> new RegionWithSequence(
                        region,
                        mRefGenome.getBases(region.chromosome(), region.start(), region.end())))
                .filter(region -> isSequenceNormal(region.sequence))
                .toList();
        int denormalWindows = batch.size() - preprocessedRegions.size();
        GU_LOGGER.debug("Skipped {} windows with denormal bases", denormalWindows);

        GU_LOGGER.debug("Running BWA-MEM alignment");
        List<byte[]> sequences = preprocessedRegions.stream().map(RegionWithSequence::sequence).toList();
        List<List<BwaMemAlignment>> alignments = mAligner.alignSeqs(sequences);
        if (alignments.size() != sequences.size()) {
            // Presumably this shouldn't occur, but we'll check to give a nicer error just in case.
            throw new RuntimeException("Alignment failed");
        }

        GU_LOGGER.debug("Running risk model");
        List<RiskModelResult> riskModelResults = alignments.stream().map(this::computeRiskModel).toList();

        GU_LOGGER.debug("Writing results");
        for (int i = 0; i < preprocessedRegions.size(); ++i) {
            writeWindowResult(preprocessedRegions.get(i).region, riskModelResults.get(i));
        }

        return new ProcessingStats(batch.size(), preprocessedRegions.size(), denormalWindows);
    }

    private record RegionWithSequence(ChrBaseRegion region, byte[] sequence) {}

    private static boolean isSequenceNormal(byte[] sequence) {
        for (byte b : sequence) {
            switch (b) {
                case 'A':
                case 'C':
                case 'G':
                case 'T':
                case 'a':
                case 'c':
                case 'g':
                case 't':
                    continue;
                default:
                    return false;
            }
        }
        return true;
    }

    private record RiskModelResult(
        // Range: [0, 1].
        // Roughly reciprocal to the number of exact match alignments. E.g. 1 means no off-target, 0.5 means 1 off-target.
        double qualityScore,
        // Range: [0, +inf].
        // Higher is higher chance of off-target match.
        long riskScore,
        // Number of alignments contributing to the risk score (i.e. above the threshold).
        int offTargetCount,
        // Raw sum of alignment scores of alignments contributing to the risk score.
        long offTargetScoreSum
    ) {}

    private RiskModelResult computeRiskModel(List<BwaMemAlignment> alignments) {
        List<Integer> offTarget = alignments.stream()
                // Only need the alignment scores. The alignment score from BWA-MEM is effectively a similarity score.
                .map(BwaMemAlignment::getAlignerScore)
                // Order by best match first.
                .sorted(Comparator.reverseOrder())
                // Drop the first alignment which is assumed to be the on-target exact match.
                .skip(1)
                // Keep only alignments with score above the configured threshold.
                .takeWhile(score -> score >= mMatchScoreThreshold)
                .toList();
        int offTargetCount = offTarget.size();
        long offTargetScoreSum = offTarget.stream().mapToLong(s -> s).sum();
        long riskScore = offTargetScoreSum - (long)offTargetCount * (mMatchScoreThreshold - mMatchScoreOffset);
        // This value is the approximate equivalent total count of exact match alignments.
        double effectiveOffTargetMatchLength = (double) riskScore / (mBaseWindowLength - mMatchScoreThreshold + mMatchScoreOffset);
        double qualityScore = 1 / (1 + effectiveOffTargetMatchLength);
        return new RiskModelResult(qualityScore, riskScore, offTargetCount, offTargetScoreSum);
    }

    private void writeWindowResult(ChrBaseRegion region, RiskModelResult risk)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(region.chromosome()).add(String.valueOf(region.start())).add(String.valueOf(region.end()))
                    .add(String.valueOf(formatQualityScore(risk.qualityScore)));
            if (mVerboseOutput) {
                sj.add(String.valueOf(risk.riskScore));
                sj.add(String.valueOf(risk.offTargetCount));
                sj.add(String.valueOf(risk.offTargetScoreSum));
            }
            mOutputWriter.write(sj.toString());
            mOutputWriter.newLine();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write output file: {}", e.toString());
        }
    }

    private String formatQualityScore(double qualityScore) {
        DecimalFormat format = new DecimalFormat();
        format.setMinimumFractionDigits(0);
        format.setMaximumFractionDigits(QUALITY_SCORE_PRECISION);
        return format.format(qualityScore);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        addRefGenomeFile(configBuilder, true);
        SpecificRegions.addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.addPath(BWA_LIB_PATH, false, BWA_LIB_PATH_DESC);

        configBuilder.addInteger(CFG_BASE_WINDOW_LENGTH, "Base window length for analysis", DEFAULT_WINDOW_LENGTH);
        configBuilder.addInteger(CFG_BASE_WINDOW_SPACING, "Offset through the genome of each base window", DEFAULT_WINDOW_SPACING);
        configBuilder.addInteger(CFG_MATCH_SCORE_THRESHOLD, "Minimum alignment score to consider a match against the window", DEFAULT_MATCH_SCORE_THRESHOLD);
        configBuilder.addInteger(CFG_MATCH_SCORE_OFFSET, "Points contributed when alignment score = threshold", DEFAULT_MATCH_SCORE_OFFSET);

        configBuilder.addInteger(CFG_BATCH_SIZE, "Number of windows to align simultaneously", DEFAULT_BATCH_SIZE);
        TaskExecutor.addThreadOptions(configBuilder);

        configBuilder.addConfigItem(CFG_OUTPUT_FILE, true, "Output filename");
        configBuilder.addFlag(CFG_VERBOSE_OUTPUT, "Output more risk info (useful for debugging)");

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ProbeQualityProfiler probeQualityProfiler = new ProbeQualityProfiler(configBuilder);
        probeQualityProfiler.run();
    }
}

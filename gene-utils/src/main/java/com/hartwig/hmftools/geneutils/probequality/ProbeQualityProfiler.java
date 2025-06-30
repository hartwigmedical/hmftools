package com.hartwig.hmftools.geneutils.probequality;

import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH_DESC;
import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
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
import java.util.List;
import java.util.StringJoiner;
import java.util.function.Supplier;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.perf.TaskExecutor;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

// Tool for helping with panel probe design. It produces a file which annotates the genome with information informing how likely a probe is
// to hybridise with off-target genome regions.
// The genome is partitioned into small base windows which are scored according to a probe off-target risk model.
// The model uses BWA-MEM alignments as a heuristic for the chance of off-target hybridisation, based on the assumption that hybridisation
// occurs with a sufficiently long exact match. The model was developed and validated against data from a previous panel design.
// Output is a TSV file containing scores for each window in the genome.
public class ProbeQualityProfiler
{
    private final BaseWindowGenerator mBaseWindowGenerator;

    private final ProbeQualityModel mProbeQualityModel;

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
        final RefGenomeSource mRefGenome = loadRefGenome(refGenome);

        SpecificRegions specificRegions = SpecificRegions.from(configBuilder);
        if (specificRegions == null) {
            throw new RuntimeException("Invalid config");
        }

        int baseWindowLength = configBuilder.getInteger(CFG_BASE_WINDOW_LENGTH);
        if (baseWindowLength < 15) {
            // Less than 15 bases probably doesn't make much sense and will cause issues trying to find appropriate params for BWA-MEM.
            throw new RuntimeException(String.format("%s must be >= 10", CFG_BASE_WINDOW_LENGTH));
        }
        GU_LOGGER.debug("Base window length: {}", baseWindowLength);

        int baseWindowSpacing = configBuilder.getInteger(CFG_BASE_WINDOW_SPACING);
        if (baseWindowSpacing < 1) {
            throw new RuntimeException(String.format("%s must be >= 1", CFG_BASE_WINDOW_SPACING));
        }
        GU_LOGGER.debug("Base window spacing: {}", baseWindowSpacing);

        int batchSize = configBuilder.getInteger(CFG_BATCH_SIZE);
        if (batchSize < 1) {
            throw new RuntimeException(String.format("%s must be >= 1", CFG_BATCH_SIZE));
        }
        GU_LOGGER.debug("Batch size: {}", batchSize);

        int matchScoreThreshold = configBuilder.getInteger(CFG_MATCH_SCORE_THRESHOLD);
        if (matchScoreThreshold > baseWindowLength) {
            // If this is true then all alignments will be excluded which is useless.
            throw new RuntimeException(String.format("%s must be <= %s", CFG_MATCH_SCORE_THRESHOLD, CFG_BASE_WINDOW_LENGTH));
        }
        GU_LOGGER.debug("Match score threshold: {}", matchScoreThreshold);

        int matchScoreOffset = configBuilder.getInteger(CFG_MATCH_SCORE_OFFSET);
        GU_LOGGER.debug("Match score offset: {}", matchScoreOffset);

        int threads = TaskExecutor.parseThreads(configBuilder);
        if (threads < 1) {
            throw new RuntimeException(String.format("%s must be >= 1", TaskExecutor.THREADS));
        }
        GU_LOGGER.debug("Threads: {}", threads);

        String refGenomeImageFile = refGenome + ".img";
        loadAlignerLibrary(configBuilder.getValue(BWA_LIB_PATH));
        Supplier<BwaMemAligner> alignerFactory = () -> createAligner(refGenomeImageFile, threads);

        mVerboseOutput = configBuilder.hasFlag(CFG_VERBOSE_OUTPUT);
        GU_LOGGER.debug("Verbose output: {}", mVerboseOutput);

        String outputFile = configBuilder.getValue(CFG_OUTPUT_FILE);
        GU_LOGGER.debug("Output file: {}", outputFile);

        mBaseWindowGenerator = new BaseWindowGenerator(mRefGenome, specificRegions, baseWindowLength, baseWindowSpacing, batchSize);
        mProbeQualityModel = new ProbeQualityModel(alignerFactory, baseWindowLength, matchScoreThreshold, matchScoreOffset);
        mOutputWriter = initialiseOutputWriter(outputFile, mVerboseOutput);
    }

    private static BwaMemAligner createAligner(String refGenomeImageFile, int threads)
    {
        if(!Files.exists(Paths.get(refGenomeImageFile)) || refGenomeImageFile.isEmpty())
        {
            throw new RuntimeException("Reference genome file is missing or empty");
        }

        BwaMemIndex index = new BwaMemIndex(refGenomeImageFile);
        BwaMemAligner aligner = new BwaMemAligner(index);

        aligner.setNThreadsOption(threads);

        return aligner;
    }

    private static BufferedWriter initialiseOutputWriter(String path, boolean verboseOutput)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(path, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME).add(FLD_POSITION_START).add(FLD_POSITION_END)
                    .add(FLD_QUALITY_SCORE);
            if (verboseOutput) {
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

        ProcessingStats stats = processBaseWindowBatches(mBaseWindowGenerator.createBaseWindowBatches());

        closeBufferedWriter(mOutputWriter);

        GU_LOGGER.info("Analysis complete, mins({})", runTimeMinsStr(startTimeMs));
        GU_LOGGER.info("Window stats:");
        GU_LOGGER.info("  Total: {}", stats.totalWindows);
        GU_LOGGER.info("  Analysed: {}", stats.totalWindows);
        GU_LOGGER.info("  Denormal: {}", stats.denormalWindows);
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

    private ProcessingStats processBaseWindowBatches(Stream<List<BaseWindowGenerator.BaseWindow>> baseWindows) {
        GU_LOGGER.info("Processing base windows");
        return baseWindows
                .map(this::processBaseWindowBatch)
                .reduce(new ProcessingStats(), ProcessingStats::add);
    }

    private ProcessingStats processBaseWindowBatch(List<BaseWindowGenerator.BaseWindow> batch) {
        GU_LOGGER.debug("Processing base window batch of size {}", batch.size());
        GU_LOGGER.debug("First base window: {}", batch.get(0).region().toString());

        GU_LOGGER.debug("Retrieving base window sequences");
        // Somewhat awkward stream processing here because we need the base sequence to filter on
        // and also need to keep the region info for later.
        List<BaseWindowGenerator.BaseWindow> filteredWindows = batch.stream()
                .filter(window -> BaseWindowGenerator.isSequenceNormal(window.sequence()))
                .toList();
        int denormalWindows = batch.size() - filteredWindows.size();
        GU_LOGGER.debug("Skipped {} windows with denormal bases", denormalWindows);

        List<byte[]> sequences = filteredWindows.stream().map(BaseWindowGenerator.BaseWindow::sequence).toList();
        List<ProbeQualityModel.Result> modelResults = mProbeQualityModel.compute(sequences);

        GU_LOGGER.debug("Writing results");
        for (int i = 0; i < filteredWindows.size(); ++i) {
            writeBaseWindowResult(filteredWindows.get(i).region(), modelResults.get(i));
        }

        return new ProcessingStats(batch.size(), filteredWindows.size(), denormalWindows);
    }

    private void writeBaseWindowResult(ChrBaseRegion region, ProbeQualityModel.Result result)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(region.chromosome()).add(String.valueOf(region.start())).add(String.valueOf(region.end()))
                    .add(String.valueOf(formatQualityScore(result.qualityScore())));
            if (mVerboseOutput) {
                sj.add(String.valueOf(result.riskScore()));
                sj.add(String.valueOf(result.offTargetCount()));
                sj.add(String.valueOf(result.offTargetScoreSum()));
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
        configBuilder.addInteger(CFG_MATCH_SCORE_OFFSET, "Risk points contributed when alignment score = threshold", DEFAULT_MATCH_SCORE_OFFSET);

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

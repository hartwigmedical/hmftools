package com.hartwig.hmftools.panelbuilder.probequality;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH_DESC;
import static com.hartwig.hmftools.common.bwa.BwaUtils.LIBBWA_PATH;
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
import static com.hartwig.hmftools.panelbuilder.probequality.Utils.createBwaMemAligner;

import java.io.BufferedWriter;
import java.io.IOException;
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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;

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
    // If true, output more information than usual. Useful for debugging.
    private final boolean mVerboseOutput;
    // Number of decimals included in the quality score.
    private static final int QUALITY_SCORE_PRECISION = 2;
    private static final int QUALITY_SCORE_PRECISION_VERBOSE = 6;

    private static final String BWA_INDEX_IMAGE_FILE_CONFIG = "bwa_index_image";
    private static final String BWA_INDEX_IMAGE_FILE_DESC = "Reference genome BWA-MEM index GATK image file";

    private static final String BASE_WINDOW_LENGTH_CONFIG = "window_length";
    private static final String BASE_WINDOW_LENGTH_DESC = "Base window length for analysis";
    private static final int BASE_WINDOW_LENGTH_DEFAULT = 40;

    private static final String BASE_WINDOW_SPACING_CONFIG = "window_spacing";
    private static final String BASE_WINDOW_SPACING_DESC = "Offset through the genome of each base window";
    private static final int BASE_WINDOW_SPACING_DEFAULT = 20;

    private static final String MATCH_SCORE_THRESHOLD_CONFIG = "match_score_threshold";
    private static final String MATCH_SCORE_THRESHOLD_DESC = "Minimum alignment score to consider a match against the window";
    private static final int MATCH_SCORE_THRESHOLD_DEFAULT = 24;    // Determined empirically via experiment

    private static final String MATCH_SCORE_OFFSET_CONFIG = "match_score_offset";
    private static final String MATCH_SCORE_OFFSET_DESC = "Risk points contributed when alignment score = threshold";
    private static final int MATCH_SCORE_OFFSET_DEFAULT = 26;   // Determined empirically via experiment

    private static final String BATCH_SIZE_CONFIG = "batch_size";
    private static final String BATCH_SIZE_DESC = "Number of windows to align simultaneously";
    private static final int BATCH_SIZE_DEFAULT = 25000;

    private static final String OUTPUT_FILE_CONFIG = "output_file";
    private static final String OUTPUT_FILE_DESC = "Output filename";

    private static final String VERBOSE_OUTPUT_CONFIG = "verbose_output";
    private static final String VERBOSE_OUTPUT_DESC = "Output more risk info (useful for debugging)";

    private static final String QUALITY_SCORE_FIELD = "QualityScore";
    private static final String RISK_SCORE_FIELD = "RiskScore";
    private static final String OFF_TARGET_COUNT_FIELD = "OffTargetCount";
    private static final String OFF_TARGET_SCORE_SUM_FIELD = "OffTargetScoreSum";

    private static final Logger LOGGER = LogManager.getLogger(ProbeQualityProfiler.class);

    public ProbeQualityProfiler(final ConfigBuilder configBuilder)
    {
        String refGenomePath = configBuilder.getValue(REF_GENOME);
        LOGGER.debug("Ref genome: {}", refGenomePath);
        RefGenomeSource refGenome = loadRefGenome(refGenomePath);
        if(refGenome == null)
        {
            throw new RuntimeException("Failed to load reference genome");
        }

        SpecificRegions specificRegions = SpecificRegions.from(configBuilder);
        if(specificRegions == null)
        {
            throw new RuntimeException("Invalid config");
        }

        int baseWindowLength = configBuilder.getInteger(BASE_WINDOW_LENGTH_CONFIG);
        if(baseWindowLength < 15)
        {
            // Less than 15 bases probably doesn't make much sense and will cause issues trying to find appropriate params for BWA-MEM.
            throw new RuntimeException(format("%s must be >= 10", BASE_WINDOW_LENGTH_CONFIG));
        }
        LOGGER.debug("Base window length: {}", baseWindowLength);

        int baseWindowSpacing = configBuilder.getInteger(BASE_WINDOW_SPACING_CONFIG);
        if(baseWindowSpacing < 1)
        {
            throw new RuntimeException(format("%s must be >= 1", BASE_WINDOW_SPACING_CONFIG));
        }
        LOGGER.debug("Base window spacing: {}", baseWindowSpacing);

        int batchSize = configBuilder.getInteger(BATCH_SIZE_CONFIG);
        if(batchSize < 1)
        {
            throw new RuntimeException(format("%s must be >= 1", BATCH_SIZE_CONFIG));
        }
        LOGGER.debug("Batch size: {}", batchSize);

        int matchScoreThreshold = configBuilder.getInteger(MATCH_SCORE_THRESHOLD_CONFIG);
        if(matchScoreThreshold > baseWindowLength)
        {
            // If this is true then all alignments will be excluded which is useless.
            throw new RuntimeException(format("%s must be <= %s", MATCH_SCORE_THRESHOLD_CONFIG, BASE_WINDOW_LENGTH_CONFIG));
        }
        LOGGER.debug("Match score threshold: {}", matchScoreThreshold);

        int matchScoreOffset = configBuilder.getInteger(MATCH_SCORE_OFFSET_CONFIG);
        LOGGER.debug("Match score offset: {}", matchScoreOffset);

        int threads = TaskExecutor.parseThreads(configBuilder);
        if(threads < 1)
        {
            throw new RuntimeException(format("%s must be >= 1", TaskExecutor.THREADS));
        }
        LOGGER.debug("Threads: {}", threads);

        loadAlignerLibrary(configBuilder.getValue(BWA_LIB_PATH));
        LOGGER.debug("BWA-MEM library path: {}", System.getProperty(LIBBWA_PATH));

        String bwaIndexImageFile = configBuilder.getValue(BWA_INDEX_IMAGE_FILE_CONFIG, refGenomePath + ".img");
        LOGGER.debug("BWA-MEM index image: {}", bwaIndexImageFile);
        Supplier<BwaMemAligner> alignerFactory = () -> createBwaMemAligner(bwaIndexImageFile, threads);

        mVerboseOutput = configBuilder.hasFlag(VERBOSE_OUTPUT_CONFIG);
        LOGGER.debug("Verbose output: {}", mVerboseOutput);

        String outputFile = configBuilder.getValue(OUTPUT_FILE_CONFIG);
        LOGGER.debug("Output file: {}", outputFile);

        mBaseWindowGenerator = new BaseWindowGenerator(refGenome, specificRegions, baseWindowLength, baseWindowSpacing, batchSize);
        mProbeQualityModel = new ProbeQualityModel(alignerFactory, baseWindowLength, matchScoreThreshold, matchScoreOffset);
        mOutputWriter = initialiseOutputWriter(outputFile, mVerboseOutput);
    }

    private static BufferedWriter initialiseOutputWriter(String path, boolean verboseOutput)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(path, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME).add(FLD_POSITION_START).add(QUALITY_SCORE_FIELD);
            if(verboseOutput)
            {
                sj.add(FLD_POSITION_END);
                sj.add(RISK_SCORE_FIELD);
                sj.add(OFF_TARGET_COUNT_FIELD);
                sj.add(OFF_TARGET_SCORE_SUM_FIELD);
            }
            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    public void run()
    {
        LOGGER.info("Starting");

        long startTimeMs = System.currentTimeMillis();

        ProcessingStats stats = processBaseWindowBatches(mBaseWindowGenerator.createBaseWindowBatches());

        closeBufferedWriter(mOutputWriter);

        LOGGER.info("Analysis complete, mins({})", runTimeMinsStr(startTimeMs));
        LOGGER.info("Window stats:");
        LOGGER.info("  Total: {}", stats.totalWindows);
        LOGGER.info("  Analysed: {}", stats.totalWindows);
        LOGGER.info("  Denormal: {}", stats.denormalWindows);
    }

    private record ProcessingStats(
            // Total number of windows considered.
            long totalWindows,
            // Number of windows analysed (i.e. not skipped).
            long analysedWindows,
            // How many base windows were skipped because they contained bases which are not ACGT.
            long denormalWindows
    )
    {
        ProcessingStats()
        {
            this(0, 0, 0);
        }

        ProcessingStats add(ProcessingStats other)
        {
            return new ProcessingStats(
                    totalWindows + other.totalWindows,
                    analysedWindows + other.analysedWindows,
                    denormalWindows + other.denormalWindows
            );
        }
    }

    private ProcessingStats processBaseWindowBatches(Stream<List<BaseWindowGenerator.BaseWindow>> baseWindows)
    {
        LOGGER.info("Processing base windows");
        return baseWindows
                .map(this::processBaseWindowBatch)
                .reduce(new ProcessingStats(), ProcessingStats::add);
    }

    private ProcessingStats processBaseWindowBatch(List<BaseWindowGenerator.BaseWindow> batch)
    {
        LOGGER.debug("Processing base window batch of size {}", batch.size());
        LOGGER.debug("First base window: {}", batch.get(0).region().toString());

        LOGGER.debug("Retrieving base window sequences");
        // Somewhat awkward stream processing here because we need the base sequence to filter on
        // and also need to keep the region info for later.
        List<BaseWindowGenerator.BaseWindow> filteredWindows = batch.stream()
                .filter(window -> BaseWindowGenerator.isSequenceNormal(window.sequence()))
                .toList();
        int denormalWindows = batch.size() - filteredWindows.size();
        LOGGER.debug("Skipped {} windows with denormal bases", denormalWindows);

        List<byte[]> sequences = filteredWindows.stream().map(BaseWindowGenerator.BaseWindow::sequence).toList();
        List<ProbeQualityModel.Result> modelResults = mProbeQualityModel.computeFromSeqBytes(sequences);

        LOGGER.debug("Writing results");
        try
        {
            for(int i = 0; i < filteredWindows.size(); ++i)
            {
                writeBaseWindowResult(filteredWindows.get(i).region(), modelResults.get(i));
            }
            mOutputWriter.flush();
        }
        catch(IOException e)
        {
            throw new RuntimeException(format("Writing to output file failed: %s", e));
        }

        return new ProcessingStats(batch.size(), filteredWindows.size(), denormalWindows);
    }

    private void writeBaseWindowResult(ChrBaseRegion region, ProbeQualityModel.Result result) throws IOException
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);
        sj.add(region.chromosome()).add(String.valueOf(region.start())).add(formatQualityScore(result.qualityScore()));
        if(mVerboseOutput)
        {
            sj.add(String.valueOf(region.end()));
            sj.add(String.valueOf(result.riskScore()));
            sj.add(String.valueOf(result.offTargetCount()));
            sj.add(String.valueOf(result.offTargetScoreSum()));
        }
        mOutputWriter.write(sj.toString());
        mOutputWriter.newLine();
    }

    private String formatQualityScore(double qualityScore)
    {
        DecimalFormat format = new DecimalFormat();
        format.setMinimumFractionDigits(0);
        format.setMaximumFractionDigits(mVerboseOutput ? QUALITY_SCORE_PRECISION_VERBOSE : QUALITY_SCORE_PRECISION);
        return format.format(qualityScore);
    }

    public static void main(String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        addRefGenomeFile(configBuilder, true);
        SpecificRegions.addSpecificChromosomesRegionsConfig(configBuilder);

        configBuilder.addPath(BWA_INDEX_IMAGE_FILE_CONFIG, false, BWA_INDEX_IMAGE_FILE_DESC);

        configBuilder.addPath(BWA_LIB_PATH, false, BWA_LIB_PATH_DESC);

        configBuilder.addInteger(BASE_WINDOW_LENGTH_CONFIG, BASE_WINDOW_LENGTH_DESC, BASE_WINDOW_LENGTH_DEFAULT);
        configBuilder.addInteger(BASE_WINDOW_SPACING_CONFIG, BASE_WINDOW_SPACING_DESC, BASE_WINDOW_SPACING_DEFAULT);
        configBuilder.addInteger(MATCH_SCORE_THRESHOLD_CONFIG, MATCH_SCORE_THRESHOLD_DESC, MATCH_SCORE_THRESHOLD_DEFAULT);
        configBuilder.addInteger(MATCH_SCORE_OFFSET_CONFIG, MATCH_SCORE_OFFSET_DESC, MATCH_SCORE_OFFSET_DEFAULT);

        configBuilder.addInteger(BATCH_SIZE_CONFIG, BATCH_SIZE_DESC, BATCH_SIZE_DEFAULT);
        TaskExecutor.addThreadOptions(configBuilder);

        configBuilder.addConfigItem(OUTPUT_FILE_CONFIG, true, OUTPUT_FILE_DESC);
        configBuilder.addFlag(VERBOSE_OUTPUT_CONFIG, VERBOSE_OUTPUT_DESC);

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ProbeQualityProfiler probeQualityProfiler = new ProbeQualityProfiler(configBuilder);
        probeQualityProfiler.run();
    }
}

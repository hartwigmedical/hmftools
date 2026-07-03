package com.hartwig.hmftools.panelbuilder;

import static java.lang.System.exit;

import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH;
import static com.hartwig.hmftools.common.bwa.BwaUtils.BWA_LIB_PATH_DESC;
import static com.hartwig.hmftools.common.bwa.BwaUtils.loadAlignerLibrary;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;
import static com.hartwig.hmftools.common.perf.TaskExecutor.addThreadOptions;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputId;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkCreateOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_LENGTH;
import static com.hartwig.hmftools.panelbuilder.probequality.Utils.createBwaMemAligner;

import java.io.File;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;
import com.hartwig.hmftools.panelbuilder.probequality.ProbeQualityModel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

// Computes the PanelBuilder probe quality score for each sequence in a FASTA file, using the same ProbeQualityModel as PanelBuilder.
// Outputs a TSV of sequence label, sequence, and quality score.
public class FastaProbeQsCalc
{
    private final Config mConfig;
    private final ProbeQualityModel mProbeQualityModel;

    private static final String OUTPUT_FILE_NAME = "probe_quality.tsv";

    private static final Logger LOGGER = LogManager.getLogger(FastaProbeQsCalc.class);

    // Input FASTA sequence together with its label (the FASTA record name).
    private record LabelledSequence(String label, String sequence)
    {
    }

    private enum OutputColumns
    {
        Label,
        Sequence,
        QualityScore
    }

    public record Config(
            String inputFasta,
            String bwaIndexImageFile,
            @Nullable String bwaLibPath,
            int probeLength,
            int matchScoreThreshold,
            int matchScoreOffset,
            int threads,
            @Nullable String outputId,
            String outputDir
    )
    {
        private static final String CFG_INPUT_FASTA = "input_fasta";
        private static final String DESC_INPUT_FASTA = "Input FASTA file of sequences to score";
        private static final String CFG_BWA_INDEX_IMAGE_FILE = "bwa_index_image";
        private static final String DESC_BWA_INDEX_IMAGE_FILE = "Reference genome BWA-MEM index GATK image file";
        private static final String CFG_PROBE_LENGTH = "probe_length";
        private static final String DESC_PROBE_LENGTH = "Length of the probe sequences to score";
        private static final String CFG_MATCH_SCORE_THRESHOLD = "match_score_threshold";
        private static final String DESC_MATCH_SCORE_THRESHOLD = "Minimum alignment score to consider an off-target match";
        private static final String CFG_MATCH_SCORE_OFFSET = "match_score_offset";
        private static final String DESC_MATCH_SCORE_OFFSET = "Risk points contributed when alignment score = threshold";

        // Defaults matching the values used to build the PanelBuilder probe quality profile, so scores are consistent by default.
        // Match score threshold and offset were determined empirically via experiment.
        private static final int MATCH_SCORE_THRESHOLD_DEFAULT = 24;
        private static final int MATCH_SCORE_OFFSET_DEFAULT = 26;

        public static Config fromConfigBuilder(final ConfigBuilder configBuilder)
        {
            return new Config(
                    configBuilder.getValue(CFG_INPUT_FASTA),
                    configBuilder.getValue(CFG_BWA_INDEX_IMAGE_FILE),
                    configBuilder.getValue(BWA_LIB_PATH),
                    configBuilder.getInteger(CFG_PROBE_LENGTH),
                    configBuilder.getInteger(CFG_MATCH_SCORE_THRESHOLD),
                    configBuilder.getInteger(CFG_MATCH_SCORE_OFFSET),
                    parseThreads(configBuilder),
                    configBuilder.getValue(OUTPUT_ID),
                    parseOutputDir(configBuilder));
        }

        public static void registerConfig(final ConfigBuilder configBuilder)
        {
            configBuilder.addPath(CFG_INPUT_FASTA, true, DESC_INPUT_FASTA);
            configBuilder.addPath(CFG_BWA_INDEX_IMAGE_FILE, true, DESC_BWA_INDEX_IMAGE_FILE);
            configBuilder.addPath(BWA_LIB_PATH, false, BWA_LIB_PATH_DESC);
            configBuilder.addInteger(CFG_PROBE_LENGTH, DESC_PROBE_LENGTH, PROBE_LENGTH);
            configBuilder.addInteger(CFG_MATCH_SCORE_THRESHOLD, DESC_MATCH_SCORE_THRESHOLD, MATCH_SCORE_THRESHOLD_DEFAULT);
            configBuilder.addInteger(CFG_MATCH_SCORE_OFFSET, DESC_MATCH_SCORE_OFFSET, MATCH_SCORE_OFFSET_DEFAULT);
            addThreadOptions(configBuilder);
            configBuilder.addConfigItem(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
            addOutputId(configBuilder);
            addLoggingOptions(configBuilder);
        }
    }

    public FastaProbeQsCalc(final Config config)
    {
        mConfig = config;

        LOGGER.info("Loading prerequisite data");

        loadAlignerLibrary(mConfig.bwaLibPath());
        Supplier<BwaMemAligner> alignerFactory = () -> createBwaMemAligner(mConfig.bwaIndexImageFile(), mConfig.threads());
        mProbeQualityModel = new ProbeQualityModel(
                alignerFactory, mConfig.probeLength(), mConfig.matchScoreThreshold(), mConfig.matchScoreOffset());
    }

    public void run()
    {
        LOGGER.info("Starting FastaProbeQsCalc");

        LOGGER.debug("Config: {}", mConfig);

        long startTimeMs = System.currentTimeMillis();

        List<LabelledSequence> sequences = loadSequences(mConfig.inputFasta());
        LOGGER.info("Loaded {} sequences", sequences.size());
        validateSequenceLengths(sequences);

        LOGGER.info("Computing probe quality scores");
        List<ProbeQualityModel.Result> results =
                mProbeQualityModel.computeFromSeqString(sequences.stream().map(LabelledSequence::sequence).toList());

        checkCreateOutputDir(mConfig.outputDir());
        writeOutput(sequences, results);

        LOGGER.info("FastaProbeQsCalc complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private static List<LabelledSequence> loadSequences(final String fastaFile)
    {
        List<LabelledSequence> sequences = new ArrayList<>();
        try(FastaSequenceFile reader = new FastaSequenceFile(new File(fastaFile), true))
        {
            ReferenceSequence record;
            while((record = reader.nextSequence()) != null)
            {
                sequences.add(new LabelledSequence(record.getName(), record.getBaseString()));
            }
        }
        return sequences;
    }

    // The probe quality model requires every sequence to have the configured probe length, so reject the input early with a clear message.
    private void validateSequenceLengths(final List<LabelledSequence> sequences)
    {
        List<String> invalidLabels = sequences.stream()
                .filter(sequence -> sequence.sequence().length() != mConfig.probeLength())
                .map(LabelledSequence::label)
                .toList();
        if(!invalidLabels.isEmpty())
        {
            throw new UserInputError(String.format(
                    "%d sequence(s) are not the required length of %d bases: %s",
                    invalidLabels.size(), mConfig.probeLength(), invalidLabels));
        }
    }

    private void writeOutput(final List<LabelledSequence> sequences, final List<ProbeQualityModel.Result> results)
    {
        String outputFile = outputFilePath();
        LOGGER.info("Writing {} quality scores to {}", results.size(), outputFile);
        try(DelimFileWriter<Integer> writer = new DelimFileWriter<>(
                outputFile, OutputColumns.values(), (index, row) ->
        {
            LabelledSequence sequence = sequences.get(index);
            row.set(OutputColumns.Label, sequence.label());
            row.set(OutputColumns.Sequence, sequence.sequence());
            row.set(OutputColumns.QualityScore, results.get(index).qualityScore());
        }))
        {
            for(int i = 0; i < sequences.size(); i++)
            {
                writer.writeRow(i);
            }
        }
    }

    private String outputFilePath()
    {
        String fileName = mConfig.outputId() == null ? OUTPUT_FILE_NAME : mConfig.outputId() + "." + OUTPUT_FILE_NAME;
        return Paths.get(mConfig.outputDir(), fileName).toString();
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();

        Config.registerConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        try
        {
            Config config = Config.fromConfigBuilder(configBuilder);
            FastaProbeQsCalc application = new FastaProbeQsCalc(config);
            application.run();
        }
        catch(UserInputError e)
        {
            LOGGER.error("Bad input data: {}", e.getMessage());
            exit(1);
        }
        catch(RuntimeException e)
        {
            LOGGER.error("Runtime error", e);
            exit(1);
        }
    }
}

package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEPTIDE_LENGTHS;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class BinderConfig
{
    public final String ValidationDataFile;

    public final String ScoreFileDir;
    public final String ScoreFileId;

    // the following reference scoring files can be read using the reference files id instead of specified individually
    public final String TrainingDataFile;
    public final String PosWeightsFile; // file with computed and cached binding matrix per allele
    public final String BindLikelihoodFile;
    public final String FlankPosWeightsFile;

    public final CalcConstants Constants;
    public final boolean CalcPairs;
    public final boolean ApplyFlanks;

    public final RandomPeptideConfig RandomPeptides;

    public final String OutputDir;
    public final String OutputId;
    public final boolean WriteFrequencyData;
    public final boolean WritePosWeightMatrix;
    public final boolean WriteBindCounts;
    public final PeptideWriteType WritePeptideType;
    public final boolean WriteSummaryData;
    public final boolean WriteLikelihood;
    public final boolean WritePanLengthDistribution;
    public final boolean RunValidation;

    public final List<String> LogCalcAlleles;

    // the set of alleles which ensured to exist in the training output data (ie PWM values)
    // including for alleles which don't exist in the training data but which are required to score novel alleles
    public final List<String> RequiredOutputAlleles;

    // set of peptide lengths processed from training data and ensured to exist in training output data
    public final List<Integer> RequiredPeptideLengths;

    public static final String SCORE_FILE_ID = "score_file_id";
    public static final String SCORE_FILE_DIR = "score_file_dir";

    private static final String TRAINING_DATA_FILE = "training_data_file";
    private static final String VALIDATION_DATA_FILE = "validation_data_file";

    public static final String FILE_ID_POS_WEIGHT = "pos_weight";
    public static final String FILE_ID_LIKELIHOOD = "likelihood";
    public static final String FILE_ID_RAND_DIST = "rand_dist";
    public static final String FILE_ID_LIKELIHOOD_RAND_DIST = "likelihood_rand_dist";
    public static final String FILE_ID_FLANK_POS_WEIGHT = "flank_pos_weight";

    private static final String APPLY_FLANKS = "apply_flanks";

    private static final String WRITE_PW_MATRIX = "write_pw_matrix";
    private static final String WRITE_BIND_COUNTS = "write_bind_counts";
    private static final String WRITE_FREQ_DATA = "write_freq_data";
    private static final String WRITE_PEPTIDE_TYPE = "write_peptide_type";
    private static final String WRITE_SUMMARY_DATA = "write_summary_data";
    private static final String WRITE_PAIRS_DATA = "write_pairs";
    private static final String WRITE_LIKELIHOOD = "write_likelihood";
    public static final String WRITE_PAN_LENGTH_DIST = "write_pan_length_dist";
    private static final String RUN_VALIDATION = "run_validation";
    private static final String LOG_CALC_ALLELES = "log_calc_alleles";

    public static final String REQUIRED_OUTPUT_ALLELES = "required_output_alleles";
    private static final String REQUIRED_PEPTIDE_LENGTHS = "required_peptide_lengths";
    public static final String OUTPUT_ID = "output_id";
    public static final String THREADS = "threads";

    public BinderConfig(final CommandLine cmd)
    {
        TrainingDataFile = cmd.getOptionValue(TRAINING_DATA_FILE);
        ValidationDataFile = cmd.getOptionValue(VALIDATION_DATA_FILE);

        ScoreFileDir = cmd.hasOption(SCORE_FILE_DIR) ? checkAddDirSeparator(cmd.getOptionValue(SCORE_FILE_DIR)) : null;
        ScoreFileId = cmd.getOptionValue(SCORE_FILE_ID);

        // load reference files either by specific name or using the scoring data dir and file id
        PosWeightsFile = getScoringFilename(cmd, ScoreFileDir, ScoreFileId, FILE_ID_POS_WEIGHT);
        FlankPosWeightsFile = getScoringFilename(cmd, ScoreFileDir, ScoreFileId, FILE_ID_FLANK_POS_WEIGHT);
        BindLikelihoodFile = getScoringFilename(cmd, ScoreFileDir, ScoreFileId, FILE_ID_LIKELIHOOD);

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        Constants = new CalcConstants(cmd);

        RequiredOutputAlleles = loadRequiredOutputAlleles(cmd.getOptionValue(REQUIRED_OUTPUT_ALLELES));

        RequiredPeptideLengths = Lists.newArrayList();

        if(cmd.hasOption(REQUIRED_PEPTIDE_LENGTHS))
        {
            String[] pepLenStrings = cmd.getOptionValue(REQUIRED_PEPTIDE_LENGTHS).split(ITEM_DELIM, -1);

            for(String pepLenStr : pepLenStrings)
            {
                int peptideLength = Integer.parseInt(pepLenStr);

                if(peptideLength < MIN_PEPTIDE_LENGTH || peptideLength > REF_PEPTIDE_LENGTH)
                {
                    NE_LOGGER.error("ignoring invalid configured peptide length({})", peptideLength);
                    continue;
                }

                RequiredPeptideLengths.add(peptideLength);
            }

            NE_LOGGER.info("requiring {} peptide lengths: {}", RequiredPeptideLengths.size(), RequiredPeptideLengths);
        }
        else
        {
            RequiredPeptideLengths.addAll(DEFAULT_PEPTIDE_LENGTHS);
        }

        CalcPairs = cmd.hasOption(WRITE_PAIRS_DATA);
        ApplyFlanks = cmd.hasOption(APPLY_FLANKS);
        RunValidation = cmd.hasOption(RUN_VALIDATION);

        RandomPeptides = new RandomPeptideConfig(cmd);

        WritePosWeightMatrix = cmd.hasOption(WRITE_PW_MATRIX);
        WriteBindCounts = cmd.hasOption(WRITE_BIND_COUNTS);
        WriteFrequencyData = cmd.hasOption(WRITE_FREQ_DATA);
        WriteSummaryData = cmd.hasOption(WRITE_SUMMARY_DATA);
        WriteLikelihood = cmd.hasOption(WRITE_LIKELIHOOD);
        WritePanLengthDistribution = cmd.hasOption(WRITE_PAN_LENGTH_DIST);
        WritePeptideType = PeptideWriteType.valueOf(cmd.getOptionValue(WRITE_PEPTIDE_TYPE, PeptideWriteType.NONE.toString()));

        LogCalcAlleles = Lists.newArrayList();

        if(cmd.hasOption(LOG_CALC_ALLELES))
        {
            LogCalcAlleles.addAll(Arrays.stream(cmd.getOptionValue(LOG_CALC_ALLELES).split(ITEM_DELIM)).collect(Collectors.toList()));
        }
    }

    public static String scoreFileConfig(final String fileType) { return fileType + "_file"; }

    public static String getScoringFilename(
            final CommandLine cmd, final String scoreFileDir, final String scoreFileId, final String fileType)
    {
        // formFilename
        String configStr = scoreFileConfig(fileType);
        String configValue = cmd.getOptionValue(configStr);

        if(scoreFileDir == null)
            return configValue;

        if(configValue != null)
            return scoreFileDir + configValue;
        else
            return formFilename(fileType, scoreFileDir, scoreFileId);
    }

    public String formOutputFilename(final String fileType)
    {
        return formFilename(fileType, OutputDir, OutputId);
    }

    public static String formFilename(final String fileType, final String dir, final String fileId)
    {
        if(fileId == null || fileId.isEmpty())
            return String.format("%sbind_%s.csv", dir, fileType);
        else
            return String.format("%sbind_%s_%s.csv", dir, fileId, fileType);
    }

    public static List<String> loadRequiredOutputAlleles(final String filename)
    {
        List<String> requiredAlleles = Lists.newArrayList();

        if(filename != null)
        {
            try
            {
                requiredAlleles.addAll(Files.readAllLines(Paths.get(filename)).stream()
                        .filter(x -> !x.equals("Allele")).collect(Collectors.toList()));

                NE_LOGGER.info("loaded {} required output alleles from file({})", requiredAlleles.size(), filename);
            }
            catch(IOException e)
            {
                NE_LOGGER.error("failed to load required alleles file({}): {}", filename, e.toString());
            }
        }

        return requiredAlleles;
    }

    public static void addCmdLineArgs(Options options)
    {
        RandomPeptideConfig.addCmdLineArgs(options);
        options.addOption(TRAINING_DATA_FILE, true, "Training data file");
        options.addOption(VALIDATION_DATA_FILE, true, "Validation data file");

        options.addOption(SCORE_FILE_ID, true, "Reference file id for scoring instead of specifying individual files");
        options.addOption(SCORE_FILE_DIR, true, "Reference file directory");
        options.addOption(scoreFileConfig(FILE_ID_POS_WEIGHT), true, "Binding position weights file");
        options.addOption(scoreFileConfig(FILE_ID_LIKELIHOOD), true, "Binding likelihood file");

        options.addOption(HLA_DEFINITIONS_FILE, true, "HLA allele definitions file");

        CalcConstants.addCmdLineArgs(options);

        options.addOption(APPLY_FLANKS, false, "Use flanks for scoring if present");
        options.addOption(WRITE_PAIRS_DATA, false, "Calculate amino-acid pairs and their coocurrence");
        options.addOption(WRITE_PW_MATRIX, false, "Write computed amino-acid + position matrix data");
        options.addOption(WRITE_BIND_COUNTS, false, "Write interim bind counts data");
        options.addOption(WRITE_FREQ_DATA, false, "Write amino-acid + position frequency data");
        options.addOption(WRITE_SUMMARY_DATA, false, "Write allele summary data including AUC");
        options.addOption(WRITE_LIKELIHOOD, false, "Write relative likelihood data");
        options.addOption(WRITE_PAN_LENGTH_DIST, false, "Write pan-peptide length distribution data");
        options.addOption(WRITE_PEPTIDE_TYPE, true, "Write peptide scores and ranks - filtered by TRAINING, LIKELY_INCORRECT, else ALL");

        options.addOption(REQUIRED_PEPTIDE_LENGTHS, true, "List of peptide-lengths separated by ';'");
        options.addOption(LOG_CALC_ALLELES, true, "Log verbose calcs for alleles separated by ';'");

        options.addOption(RUN_VALIDATION, false, "Run validation routines");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file id");
        options.addOption(LOG_DEBUG, false, "Log verbose");
    }
}

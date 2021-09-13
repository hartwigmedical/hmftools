package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.NeoCommon.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEPTIDE_LENGTHS;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.ExpressionLikelihood.EXP_LIKELIHOOD_FILE;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.scoreFileConfig;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class TrainConfig
{
    // the following reference scoring files can be read using the reference files id instead of specified individually
    public final String TrainingDataFile;
    public final String RecognitionDataFile;
    public final String ExpressionLikelihoodFile;

    public final CalcConstants Constants;
    public final boolean CalcPairs;
    public final boolean ApplyFlanks;

    public final RandomPeptideConfig RandomPeptides;

    public final String OutputDir;
    public final String OutputId;
    public final boolean WriteFrequencyData;
    public final boolean WritePosWeightMatrix;
    public final boolean WriteBindCounts;
    public final boolean WriteLikelihood;
    public final boolean WritePanLengthDistribution;
    public final boolean RunValidation;

    public final List<String> LogCalcAlleles;

    // the set of alleles which ensured to exist in the training output data (ie PWM values)
    // including for alleles which don't exist in the training data but which are required to score novel alleles
    public final List<String> RequiredOutputAlleles;

    // set of peptide lengths processed from training data and ensured to exist in training output data
    public final List<Integer> RequiredPeptideLengths;

    private static final String TRAINING_DATA_FILE = "training_data_file";
    public static final String RECOGNITION_DATA_FILE = "recognition_data_file";

    public static final String FILE_ID_POS_WEIGHT = "pos_weight";
    public static final String FILE_ID_LIKELIHOOD = "likelihood";
    public static final String FILE_ID_RAND_DIST = "rand_dist";
    public static final String FILE_ID_LIKELIHOOD_RAND_DIST = "likelihood_rand_dist";
    public static final String FILE_ID_FLANK_POS_WEIGHT = "flank_pos_weight";

    private static final String APPLY_FLANKS = "apply_flanks";

    private static final String WRITE_PW_MATRIX = "write_pw_matrix";
    private static final String WRITE_BIND_COUNTS = "write_bind_counts";
    private static final String WRITE_FREQ_DATA = "write_freq_data";
    private static final String WRITE_PAIRS_DATA = "write_pairs";
    private static final String WRITE_LIKELIHOOD = "write_likelihood";
    public static final String WRITE_PAN_LENGTH_DIST = "write_pan_length_dist";
    private static final String RUN_VALIDATION = "run_validation";
    private static final String LOG_CALC_ALLELES = "log_calc_alleles";

    public static final String REQUIRED_OUTPUT_ALLELES = "required_output_alleles";
    public static final String REQUIRED_PEPTIDE_LENGTHS = "required_peptide_lengths";

    public TrainConfig(final CommandLine cmd)
    {
        TrainingDataFile = cmd.getOptionValue(TRAINING_DATA_FILE);
        RecognitionDataFile = cmd.getOptionValue(RECOGNITION_DATA_FILE);
        ExpressionLikelihoodFile = cmd.getOptionValue(EXP_LIKELIHOOD_FILE);

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
        WriteLikelihood = cmd.hasOption(WRITE_LIKELIHOOD);
        WritePanLengthDistribution = cmd.hasOption(WRITE_PAN_LENGTH_DIST);

        LogCalcAlleles = Lists.newArrayList();

        if(cmd.hasOption(LOG_CALC_ALLELES))
        {
            LogCalcAlleles.addAll(Arrays.stream(cmd.getOptionValue(LOG_CALC_ALLELES).split(ITEM_DELIM)).collect(Collectors.toList()));
        }
    }

    public String formTrainingFilename(final String fileType)
    {
        return formTrainingFilename(OutputDir, fileType, OutputId);
    }

    public static String formTrainingFilename(final String dir, final String fileType, final String fileId)
    {
        return BindCommon.formFilename(dir, String.format("train_%s", fileType), fileId);
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
        options.addOption(TRAINING_DATA_FILE, true, "Training data file");
        options.addOption(RECOGNITION_DATA_FILE, true, "Immunogenic recognition data file");
        options.addOption(HLA_DEFINITIONS_FILE, true, "HLA allele definitions file");
        options.addOption(EXP_LIKELIHOOD_FILE, true, "Expression likelihood file");

        RandomPeptideConfig.addCmdLineArgs(options);
        CalcConstants.addCmdLineArgs(options);

        options.addOption(APPLY_FLANKS, false, "Use flanks for scoring if present");
        options.addOption(WRITE_PAIRS_DATA, false, "Calculate amino-acid pairs and their coocurrence");
        options.addOption(WRITE_PW_MATRIX, false, "Write computed amino-acid + position matrix data");
        options.addOption(WRITE_BIND_COUNTS, false, "Write interim bind counts data");
        options.addOption(WRITE_FREQ_DATA, false, "Write amino-acid + position frequency data");
        options.addOption(WRITE_LIKELIHOOD, false, "Write relative likelihood data");
        options.addOption(WRITE_PAN_LENGTH_DIST, false, "Write pan-peptide length distribution data");

        options.addOption(REQUIRED_PEPTIDE_LENGTHS, true, "List of peptide-lengths separated by ';'");
        options.addOption(LOG_CALC_ALLELES, true, "Log verbose calcs for alleles separated by ';'");

        options.addOption(RUN_VALIDATION, false, "Run validation routines");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file id");
        options.addOption(LOG_DEBUG, false, "Log verbose");
    }
}

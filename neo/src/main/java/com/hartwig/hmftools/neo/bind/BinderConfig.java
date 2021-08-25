package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_ALLELE_MOTIF_WEIGHT;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEPTIDE_LENGTHS;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEPTIDE_LENGTH_WEIGHT;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_WEIGHT_EXPONENT;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class BinderConfig
{
    public final String TrainingDataFile;
    public final String ValidationDataFile;
    public final String BindMatrixFile; // file with computed and cached binding matrix per allele
    public final String BindLikelihoodFile;

    public final CalcConstants Constants;
    public final boolean RunScoring;
    public final boolean CalcPairs;

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

    // the set of alleles which ensured to exist in the training output data (ie PWM values)
    // including for alleles which don't exist in the training data but which are required to score novel alleles
    public final List<String> RequiredOutputAlleles;

    // set of peptide lengths processed from training data and ensured to exist in training output data
    public final List<Integer> RequiredPeptideLengths;

    private static final String TRAINING_DATA_FILE = "training_data_file";
    private static final String VALIDATION_DATA_FILE = "validation_data_file";
    private static final String BIND_MATRIX_FILE = "bind_matrix_file";
    private static final String BIND_LIKELIHOOD_FILE = "bind_likelihood_file";

    private static final String WEIGHT_EXPONENT = "weight_exponent";
    private static final String LENGTH_WEIGHT = "length_weight";
    private static final String ALLELE_MOTIF_WEIGHT = "allele_motif_weight";
    private static final String NOISE_PROB = "noise_prob";
    private static final String NOISE_WEIGHT = "noise_weight";
    private static final String GLOBAL_WEIGHT = "global_weight";

    private static final String RUN_SCORING = "run_scoring";
    private static final String WRITE_PW_MATRIX = "write_pw_matrix";
    private static final String WRITE_BIND_COUNTS = "write_bind_counts";
    private static final String WRITE_FREQ_DATA = "write_freq_data";
    private static final String WRITE_PEPTIDE_TYPE = "write_peptide_type";
    private static final String WRITE_SUMMARY_DATA = "write_summary_data";
    private static final String WRITE_PAIRS_DATA = "write_pairs";
    private static final String WRITE_LIKELIHOOD = "write_likelihood";
    public static final String WRITE_PAN_LENGTH_DIST = "write_pan_length_dist";
    private static final String RUN_VALIDATION = "run_validation";

    public static final String REQUIRED_OUTPUT_ALLELES = "required_output_alleles";
    private static final String REQUIRED_PEPTIDE_LENGTHS = "required_peptide_lengths";
    public static final String OUTPUT_ID = "output_id";

    public BinderConfig(final CommandLine cmd)
    {
        TrainingDataFile = cmd.getOptionValue(TRAINING_DATA_FILE);
        ValidationDataFile = cmd.getOptionValue(VALIDATION_DATA_FILE);
        BindMatrixFile = cmd.getOptionValue(BIND_MATRIX_FILE);
        BindLikelihoodFile = cmd.getOptionValue(BIND_LIKELIHOOD_FILE);

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        Constants = new CalcConstants(
                Double.parseDouble(cmd.getOptionValue(LENGTH_WEIGHT, String.valueOf(DEFAULT_PEPTIDE_LENGTH_WEIGHT))),
                Double.parseDouble(cmd.getOptionValue(ALLELE_MOTIF_WEIGHT, String.valueOf(DEFAULT_ALLELE_MOTIF_WEIGHT))),
                Double.parseDouble(cmd.getOptionValue(WEIGHT_EXPONENT, String.valueOf(DEFAULT_WEIGHT_EXPONENT))),
                Double.parseDouble(cmd.getOptionValue(NOISE_PROB, "0")), // String.valueOf(DEFAULT_NOISE_PROB)
                Double.parseDouble(cmd.getOptionValue(NOISE_WEIGHT, "0")), // String.valueOf(DEFAULT_NOISE_WEIGHT)
                Double.parseDouble(cmd.getOptionValue(GLOBAL_WEIGHT, "0")));

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
        RunScoring = cmd.hasOption(RUN_SCORING);
        RunValidation = cmd.hasOption(RUN_VALIDATION);

        RandomPeptides = new RandomPeptideConfig(cmd);

        WritePosWeightMatrix = cmd.hasOption(WRITE_PW_MATRIX);
        WriteBindCounts = cmd.hasOption(WRITE_BIND_COUNTS);
        WriteFrequencyData = cmd.hasOption(WRITE_FREQ_DATA);
        WriteSummaryData = cmd.hasOption(WRITE_SUMMARY_DATA);
        WriteLikelihood = cmd.hasOption(WRITE_LIKELIHOOD);
        WritePanLengthDistribution = cmd.hasOption(WRITE_PAN_LENGTH_DIST);
        WritePeptideType = PeptideWriteType.valueOf(cmd.getOptionValue(WRITE_PEPTIDE_TYPE, PeptideWriteType.NONE.toString()));
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

    public String formFilename(final String fileId)
    {
        return formFilename(fileId, OutputDir, OutputId);
    }

    public static String formFilename(final String fileId, final String outputDir, final String outputId)
    {
        if(outputId == null || outputId.isEmpty())
            return String.format("%sbind_%s.csv", outputDir, fileId);
        else
            return String.format("%sbind_%s_%s.csv", outputDir, outputId, fileId);
    }

    public static void addCmdLineArgs(Options options)
    {
        RandomPeptideConfig.addCmdLineArgs(options);
        options.addOption(TRAINING_DATA_FILE, true, "Training data file");
        options.addOption(VALIDATION_DATA_FILE, true, "Validation data file");
        options.addOption(BIND_MATRIX_FILE, true, "Binding matrix data file");
        options.addOption(BIND_LIKELIHOOD_FILE, true, "Binding relative likelihood file");
        options.addOption(HLA_DEFINITIONS_FILE, true, "HLA allele definitions file");

        options.addOption(WEIGHT_EXPONENT, true, "Weight exponent");
        options.addOption(ALLELE_MOTIF_WEIGHT, true, "Allele motif weight");
        options.addOption(LENGTH_WEIGHT, true, "Length weight");
        options.addOption(NOISE_PROB, true, "Noise target probability");
        options.addOption(NOISE_WEIGHT, true, "Noise weight");
        options.addOption(GLOBAL_WEIGHT, true, "Global counts weight");

        options.addOption(RUN_SCORING, false, "Use binding matrix data to score training and random peptide data");
        options.addOption(WRITE_PAIRS_DATA, false, "Calculate amino-acid pairs and their coocurrence");
        options.addOption(WRITE_PW_MATRIX, false, "Write computed amino-acid + position matrix data");
        options.addOption(WRITE_BIND_COUNTS, false, "Write interim bind counts data");
        options.addOption(WRITE_FREQ_DATA, false, "Write amino-acid + position frequency data");
        options.addOption(WRITE_SUMMARY_DATA, false, "Write allele summary data including AUC");
        options.addOption(WRITE_LIKELIHOOD, false, "Write relative likelihood data");
        options.addOption(WRITE_PAN_LENGTH_DIST, false, "Write pan-peptide length distribution data");
        options.addOption(WRITE_PEPTIDE_TYPE, true, "Write peptide scores and ranks - filtered by TRAINING, LIKELY_INCORRECT, else ALL");

        options.addOption(REQUIRED_PEPTIDE_LENGTHS, true, "List of peptide-lengths separated by ';'");

        options.addOption(RUN_VALIDATION, false, "Run validation routines");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file id");
        options.addOption(LOG_DEBUG, false, "Log verbose");
    }
}

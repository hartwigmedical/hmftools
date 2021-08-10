package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_ALLELE_MOTIF_WEIGHT;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_BINDING_AFFINITY_HIGH;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_BINDING_AFFINITY_LOW;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_MAX_AFFINITY;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEPTIDE_LENGTH_WEIGHT;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_WEIGHT_EXPONENT;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class BinderConfig
{
    public final String TrainingDataFile;
    public final String ValidationDataFile;
    public final String BindMatrixFile; // file with computed and cached binding matrix per allele
    public final String RandomPeptidePredictionsFile; // predictions per allele for random peptides

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

    public final List<String> RequiredAlleles;
    public final List<Integer> RequiredPeptideLengths;

    private static final String TRAINING_DATA_FILE = "training_data_file";
    private static final String VALIDATION_DATA_FILE = "validation_data_file";
    private static final String RANDOM_PEPTIDE_PRED_FILE = "random_peptide_pred_file";
    private static final String BIND_MATRIX_FILE = "bind_matrix_file";

    private static final String MAX_AFFINITY = "max_affinity";
    private static final String BINDING_AFFINITY_LOW = "binding_affinity_low";
    private static final String BINDING_AFFINITY_HIGH = "binding_aff_high";
    private static final String APPLY_SCALED_COUNT = "apply_scaled";
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

    private static final String REQUIRED_ALLELES = "required_alleles";
    private static final String REQUIRED_PEPTIDE_LENGTHS = "required_peptide_lengths";
    public static final String OUTPUT_ID = "output_id";

    public BinderConfig(final CommandLine cmd)
    {
        TrainingDataFile = cmd.getOptionValue(TRAINING_DATA_FILE);
        ValidationDataFile = cmd.getOptionValue(VALIDATION_DATA_FILE);
        RandomPeptidePredictionsFile = cmd.getOptionValue(RANDOM_PEPTIDE_PRED_FILE);
        BindMatrixFile = cmd.getOptionValue(BIND_MATRIX_FILE);

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        Constants = new CalcConstants(
                Double.parseDouble(cmd.getOptionValue(LENGTH_WEIGHT, String.valueOf(DEFAULT_PEPTIDE_LENGTH_WEIGHT))),
                Double.parseDouble(cmd.getOptionValue(ALLELE_MOTIF_WEIGHT, String.valueOf(DEFAULT_ALLELE_MOTIF_WEIGHT))),
                Double.parseDouble(cmd.getOptionValue(WEIGHT_EXPONENT, String.valueOf(DEFAULT_WEIGHT_EXPONENT))),
                Double.parseDouble(cmd.getOptionValue(NOISE_PROB, "0")), // String.valueOf(DEFAULT_NOISE_PROB)
                Double.parseDouble(cmd.getOptionValue(NOISE_WEIGHT, "0")), // String.valueOf(DEFAULT_NOISE_WEIGHT)
                Double.parseDouble(cmd.getOptionValue(GLOBAL_WEIGHT, "0")),
                Double.parseDouble(cmd.getOptionValue(MAX_AFFINITY, String.valueOf(DEFAULT_MAX_AFFINITY))),
                Double.parseDouble(cmd.getOptionValue(BINDING_AFFINITY_LOW, String.valueOf(DEFAULT_BINDING_AFFINITY_LOW))),
                Double.parseDouble(cmd.getOptionValue(BINDING_AFFINITY_HIGH, String.valueOf(DEFAULT_BINDING_AFFINITY_HIGH))),
                cmd.hasOption(APPLY_SCALED_COUNT));

        RequiredAlleles = Lists.newArrayList();

        if(cmd.hasOption(REQUIRED_ALLELES))
        {
            String requiredAlleles = cmd.getOptionValue(REQUIRED_ALLELES);

            if(requiredAlleles.endsWith(".csv") && Files.exists(Paths.get(requiredAlleles)))
            {
                try
                {
                    Files.readAllLines(Paths.get(requiredAlleles)).stream()
                            .filter(x -> !x.equals("Allele")).forEach(x -> RequiredAlleles.add(x));

                    NE_LOGGER.info("loaded {} required alleles from file({})", RequiredAlleles.size(), requiredAlleles);
                }
                catch(IOException e)
                {
                    NE_LOGGER.error("failed to load required alleles file({}): {}", requiredAlleles, e.toString());
                }
            }
            else
            {
                Arrays.stream(cmd.getOptionValue(REQUIRED_ALLELES).split(ITEM_DELIM, -1)).forEach(x -> RequiredAlleles.add(x));
                NE_LOGGER.info("filtering for {} alleles: {}", RequiredAlleles.size(), RequiredAlleles);
            }
        }

        RequiredPeptideLengths = Lists.newArrayList();

        if(cmd.hasOption(REQUIRED_PEPTIDE_LENGTHS))
        {
            Arrays.stream(cmd.getOptionValue(REQUIRED_PEPTIDE_LENGTHS).split(ITEM_DELIM, -1))
                    .forEach(x -> RequiredPeptideLengths.add(Integer.parseInt(x)));
            NE_LOGGER.info("requiring {} peptide lengths: {}", RequiredPeptideLengths.size(), RequiredPeptideLengths);
        }

        CalcPairs = cmd.hasOption(WRITE_PAIRS_DATA);
        RunScoring = cmd.hasOption(RUN_SCORING);

        RandomPeptides = new RandomPeptideConfig(cmd);

        WritePosWeightMatrix = cmd.hasOption(WRITE_PW_MATRIX);
        WriteBindCounts = cmd.hasOption(WRITE_BIND_COUNTS);
        WriteFrequencyData = cmd.hasOption(WRITE_FREQ_DATA);
        WriteSummaryData = cmd.hasOption(WRITE_SUMMARY_DATA);
        WriteLikelihood = cmd.hasOption(WRITE_LIKELIHOOD);
        WritePeptideType = PeptideWriteType.valueOf(cmd.getOptionValue(WRITE_PEPTIDE_TYPE, PeptideWriteType.NONE.toString()));
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
        options.addOption(RANDOM_PEPTIDE_PRED_FILE, true, "Random peptide predictions file");
        options.addOption(BIND_MATRIX_FILE, true, "Binding matrix data file");
        options.addOption(HLA_DEFINITIONS_FILE, true, "HLA allele definitions file");

        options.addOption(BINDING_AFFINITY_HIGH, true, "Upper binding affinity threshold");
        options.addOption(BINDING_AFFINITY_LOW, true, "Lower binding affinity threshold");
        options.addOption(MAX_AFFINITY, true, "Binding affinity exponent  for score calc: 1 - log(exp,affinity)");
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
        options.addOption(WRITE_PEPTIDE_TYPE, true, "Write peptide scores and ranks - filtered by TRAINING, LIKELY_INCORRECT, else ALL");
        options.addOption(APPLY_SCALED_COUNT, false, "Calculate amino-acid pairs and their coocurrence");

        options.addOption(REQUIRED_ALLELES, true, "List of alleles separated by ';'");
        options.addOption(REQUIRED_PEPTIDE_LENGTHS, true, "List of peptide-lengths separated by ';'");

        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file id");
        options.addOption(LOG_DEBUG, false, "Log verbose");
    }
}

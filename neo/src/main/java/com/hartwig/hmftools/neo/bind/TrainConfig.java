package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.IGNORE_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.CSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindCommon.FLD_ALLELE;
import static com.hartwig.hmftools.neo.bind.BindConstants.DEFAULT_PEPTIDE_LENGTHS;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.ExpressionLikelihood.EXP_LIKELIHOOD_FILE;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

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
    public static final String FILE_ID_FLANK_POS_WEIGHT = "flank_pos_weight";
    public static final String FILE_ID_RECOGNITION = "recognition";
    public static final String FILE_ID_EXPRESSION_DIST = "expression_dist";

    public static final String FILE_ID_RAND_DIST = "rand_dist";
    public static final String FILE_ID_LIKELIHOOD_DIST = "likelihood_rand_dist";
    public static final String FILE_ID_EXP_LIKELIHOOD_DIST = "exp_likelihood_rand_dist";

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

    public TrainConfig(final ConfigBuilder configBuilder)
    {
        TrainingDataFile = configBuilder.getValue(TRAINING_DATA_FILE);
        RecognitionDataFile = configBuilder.getValue(RECOGNITION_DATA_FILE);
        ExpressionLikelihoodFile = configBuilder.getValue(EXP_LIKELIHOOD_FILE);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        Constants = new CalcConstants(configBuilder);

        RequiredOutputAlleles = loadRequiredOutputAlleles(configBuilder.getValue(REQUIRED_OUTPUT_ALLELES));

        RequiredPeptideLengths = Lists.newArrayList();

        if(configBuilder.hasValue(REQUIRED_PEPTIDE_LENGTHS))
        {
            String[] pepLenStrings = configBuilder.getValue(REQUIRED_PEPTIDE_LENGTHS).split(ITEM_DELIM, -1);

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

        CalcPairs = configBuilder.hasFlag(WRITE_PAIRS_DATA);
        ApplyFlanks = configBuilder.hasFlag(APPLY_FLANKS);
        RunValidation = configBuilder.hasFlag(RUN_VALIDATION);

        RandomPeptides = RandomPeptideConfig.forWriting(configBuilder);

        WritePosWeightMatrix = configBuilder.hasFlag(WRITE_PW_MATRIX);
        WriteBindCounts = configBuilder.hasFlag(WRITE_BIND_COUNTS);
        WriteFrequencyData = configBuilder.hasFlag(WRITE_FREQ_DATA);
        WriteLikelihood = configBuilder.hasFlag(WRITE_LIKELIHOOD);
        WritePanLengthDistribution = configBuilder.hasFlag(WRITE_PAN_LENGTH_DIST);

        LogCalcAlleles = Lists.newArrayList();

        if(configBuilder.hasValue(LOG_CALC_ALLELES))
        {
            LogCalcAlleles.addAll(Arrays.stream(configBuilder.getValue(LOG_CALC_ALLELES).split(ITEM_DELIM)).collect(Collectors.toList()));
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
                List<String> lines = Files.readAllLines(Paths.get(filename));
                String header = lines.get(0);
                lines.remove(0);

                Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, CSV_DELIM);
                int alleleIndex = fieldsIndexMap.get(FLD_ALLELE);

                for(String line : lines)
                {
                    if(line.isEmpty() || line.startsWith(IGNORE_SAMPLE_ID))
                        continue;

                    String[] values = line.split(CSV_DELIM, -1);
                    String allele = values[alleleIndex];

                    // convert from Lilac naming convention if required
                    requiredAlleles.add(allele);
                }

                NE_LOGGER.info("loaded {} required output alleles from file({})", requiredAlleles.size(), filename);
            }
            catch(IOException e)
            {
                NE_LOGGER.error("failed to load required alleles file({}): {}", filename, e.toString());
            }
        }

        return requiredAlleles;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(TRAINING_DATA_FILE, true, "Training data file");
        configBuilder.addPath(RECOGNITION_DATA_FILE, false, "Immunogenic recognition data file");
        configBuilder.addPath(HLA_DEFINITIONS_FILE, true, "HLA allele definitions file");
        configBuilder.addPath(EXP_LIKELIHOOD_FILE, false, "Expression likelihood file");

        RandomPeptideConfig.addConfigForWriting(configBuilder);
        CalcConstants.addConfig(configBuilder);

        configBuilder.addFlag(APPLY_FLANKS, "Use flanks for scoring if present");
        configBuilder.addFlag(WRITE_PAIRS_DATA, "Calculate amino-acid pairs and their coocurrence");
        configBuilder.addFlag(WRITE_PW_MATRIX, "Write computed amino-acid + position matrix data");
        configBuilder.addFlag(WRITE_BIND_COUNTS, "Write interim bind counts data");
        configBuilder.addFlag(WRITE_FREQ_DATA, "Write amino-acid + position frequency data");
        configBuilder.addFlag(WRITE_LIKELIHOOD, "Write relative likelihood data");
        configBuilder.addFlag(WRITE_PAN_LENGTH_DIST, "Write pan-peptide length distribution data");

        configBuilder.addConfigItem(REQUIRED_PEPTIDE_LENGTHS, false, "List of peptide-lengths separated by ';'");
        configBuilder.addConfigItem(LOG_CALC_ALLELES, false, "Log verbose calcs for alleles separated by ';'");

        configBuilder.addFlag(RUN_VALIDATION, "Run validation routines");
        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
    }
}

package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_EXPRESSION_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_FLANK_POS_WEIGHT;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_LIKELIHOOD;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_POS_WEIGHT;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_RECOGNITION;
import static com.hartwig.hmftools.neo.bind.TrainConfig.formTrainingFilename;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class ScoreConfig
{
    // allele + peptide to score & rank using training data
    public final String ValidationDataFile;
    public final String RecognitionDataFile;
    public final boolean CheckSelfRecognition;

    public final String ScoreFileDir;
    public final String ScoreFileId;

    // the following reference scoring files can be read using the reference files id instead of specified individually
    public final String PosWeightsFile; // file with computed and cached binding matrix per allele
    public final String BindLikelihoodFile;
    public final String FlankPosWeightsFile;
    public final String ExpressionLikelihoodFile;

    public final RandomPeptideConfig RandomPeptides;

    public final boolean WritePeptideScores;
    public final boolean WriteSummaryData;

    public final String OutputDir;
    public final String OutputId;

    public static final String SCORE_FILE_ID = "score_file_id";
    public static final String SCORE_FILE_DIR = "score_file_dir";

    private static final String VALIDATION_DATA_FILE = "validation_data_file";
    private static final String CHECK_SELF_RECOGNITION = "check_self_recog";

    private static final String WRITE_PEPTIDE_SCORES = "write_peptide_scores";
    private static final String WRITE_SUMMARY_DATA = "write_summary_data";

    public ScoreConfig(final ConfigBuilder configBuilder)
    {
        ValidationDataFile = configBuilder.getValue(VALIDATION_DATA_FILE);

        ScoreFileDir = checkAddDirSeparator(configBuilder.getValue(SCORE_FILE_DIR));
        ScoreFileId = configBuilder.getValue(SCORE_FILE_ID);

        // load reference files either by specific name or using the scoring data dir and file id
        PosWeightsFile = getScoringFilename(configBuilder, ScoreFileDir, ScoreFileId, FILE_ID_POS_WEIGHT);
        FlankPosWeightsFile = getScoringFilename(configBuilder, ScoreFileDir, ScoreFileId, FILE_ID_FLANK_POS_WEIGHT);
        BindLikelihoodFile = getScoringFilename(configBuilder, ScoreFileDir, ScoreFileId, FILE_ID_LIKELIHOOD);
        RecognitionDataFile = getScoringFilename(configBuilder, ScoreFileDir, ScoreFileId, FILE_ID_RECOGNITION);
        ExpressionLikelihoodFile = getScoringFilename(configBuilder, ScoreFileDir, ScoreFileId, FILE_ID_EXPRESSION_DIST);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        RandomPeptides = RandomPeptideConfig.forReading(configBuilder);

        WriteSummaryData = configBuilder.hasFlag(WRITE_SUMMARY_DATA);
        WritePeptideScores = configBuilder.hasFlag(WRITE_PEPTIDE_SCORES);

        CheckSelfRecognition = configBuilder.hasFlag(CHECK_SELF_RECOGNITION)
                || (ValidationDataFile != null && ValidationDataFile.equals(RecognitionDataFile));
    }

    public static String scoreFileConfig(final String fileType) { return fileType + "_file"; }

    public static String getScoringFilename(
            final ConfigBuilder configBuilder, final String scoreFileDir, final String scoreFileId, final String fileType)
    {
        // formFilename
        String configStr = scoreFileConfig(fileType);
        String configValue = configBuilder.getValue(configStr);

        if(scoreFileDir == null)
            return configValue;

        if(configValue != null)
            return scoreFileDir + configValue;
        else
            return formTrainingFilename(scoreFileDir, fileType, scoreFileId);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        RandomPeptideConfig.addConfigForReading(configBuilder);

        configBuilder.addConfigItem(SCORE_FILE_ID, false, "Reference file id for scoring instead of specifying individual files");
        configBuilder.addPath(SCORE_FILE_DIR, false, "Reference file directory");
        configBuilder.addPath(scoreFileConfig(FILE_ID_POS_WEIGHT), false, "Binding position weights file");
        configBuilder.addPath(scoreFileConfig(FILE_ID_FLANK_POS_WEIGHT), false, "Binding flank weights file");
        configBuilder.addPath(scoreFileConfig(FILE_ID_LIKELIHOOD), false, "Binding likelihood file");
        configBuilder.addPath(scoreFileConfig(FILE_ID_EXPRESSION_DIST), false, "Expression likelihood file");
        configBuilder.addPath(scoreFileConfig(FILE_ID_RECOGNITION), false, "Immunogenic pHLA recognition data file");

        configBuilder.addPath(VALIDATION_DATA_FILE, false, "Validation data file");

        configBuilder.addFlag(WRITE_SUMMARY_DATA, "Write summary results per allele");
        configBuilder.addFlag(WRITE_PEPTIDE_SCORES, "Write score and rank data per peptide");
        configBuilder.addFlag(CHECK_SELF_RECOGNITION, "Don't allow recognition with same peptide");
    }
}

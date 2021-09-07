package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.ConfigUtils.LOG_DEBUG;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_FLANK_POS_WEIGHT;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_LIKELIHOOD;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_POS_WEIGHT;
import static com.hartwig.hmftools.neo.bind.TrainConfig.RECOGNITION_DATA_FILE;
import static com.hartwig.hmftools.neo.bind.TrainConfig.formFilename;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class ScoreConfig
{
    // allele + peptide to score & rank using training data
    public final String ValidationDataFile;
    public final String RecognitionDataFile;

    public final String ScoreFileDir;
    public final String ScoreFileId;

    // the following reference scoring files can be read using the reference files id instead of specified individually
    public final String PosWeightsFile; // file with computed and cached binding matrix per allele
    public final String BindLikelihoodFile;
    public final String FlankPosWeightsFile;

    public final RandomPeptideConfig RandomPeptides;

    public final boolean WritePeptideScores;
    public final boolean WriteSummaryData;

    public final String OutputDir;
    public final String OutputId;

    public static final String SCORE_FILE_ID = "score_file_id";
    public static final String SCORE_FILE_DIR = "score_file_dir";

    private static final String VALIDATION_DATA_FILE = "validation_data_file";

    private static final String WRITE_PEPTIDE_SCORES = "write_peptide_scores";
    private static final String WRITE_SUMMARY_DATA = "write_summary_data";

    public ScoreConfig(final CommandLine cmd)
    {
        ValidationDataFile = cmd.getOptionValue(VALIDATION_DATA_FILE);
        RecognitionDataFile = cmd.getOptionValue(RECOGNITION_DATA_FILE);

        ScoreFileDir = cmd.hasOption(SCORE_FILE_DIR) ? checkAddDirSeparator(cmd.getOptionValue(SCORE_FILE_DIR)) : null;
        ScoreFileId = cmd.getOptionValue(SCORE_FILE_ID);

        // load reference files either by specific name or using the scoring data dir and file id
        PosWeightsFile = getScoringFilename(cmd, ScoreFileDir, ScoreFileId, FILE_ID_POS_WEIGHT);
        FlankPosWeightsFile = getScoringFilename(cmd, ScoreFileDir, ScoreFileId, FILE_ID_FLANK_POS_WEIGHT);
        BindLikelihoodFile = getScoringFilename(cmd, ScoreFileDir, ScoreFileId, FILE_ID_LIKELIHOOD);

        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID);

        RandomPeptides = new RandomPeptideConfig(cmd);

        WriteSummaryData = cmd.hasOption(WRITE_SUMMARY_DATA);
        WritePeptideScores = cmd.hasOption(WRITE_PEPTIDE_SCORES);
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

    public static void addCmdLineArgs(Options options)
    {
        RandomPeptideConfig.addCmdLineArgs(options);
        options.addOption(VALIDATION_DATA_FILE, true, "Validation data file");
        options.addOption(RECOGNITION_DATA_FILE, true, "Immunogenic recognition data file");

        options.addOption(SCORE_FILE_ID, true, "Reference file id for scoring instead of specifying individual files");
        options.addOption(SCORE_FILE_DIR, true, "Reference file directory");
        options.addOption(scoreFileConfig(FILE_ID_POS_WEIGHT), true, "Binding position weights file");
        options.addOption(scoreFileConfig(FILE_ID_LIKELIHOOD), true, "Binding likelihood file");

        options.addOption(WRITE_SUMMARY_DATA, false, "Write summary results per allele");
        options.addOption(WRITE_PEPTIDE_SCORES, false, "Write score and rank data per peptide");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        options.addOption(OUTPUT_ID, true, "Output file id");
        options.addOption(LOG_DEBUG, false, "Log verbose");
    }
}

package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.SCORE_FILE_DIR;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.SCORE_FILE_ID;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.getScoringFilename;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.scoreFileConfig;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_EXP_LIKELIHOOD_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_LIKELIHOOD_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_RAND_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.REQUIRED_OUTPUT_ALLELES;
import static com.hartwig.hmftools.neo.bind.TrainConfig.loadRequiredOutputAlleles;

import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class RandomPeptideConfig
{
    public final String RandomPeptidesFile; // list of random peptides from the proteome
    public final String ScoreDistributionFile; // internally scored and ranked distribuion of random peptides
    public final String LikelihoodDistributionFile;
    public final String ExpressionLikelihoodDistributionFile;
    public final String OutputDir;
    public final String OutputId;

    public final boolean WriteRandomDistribution;
    public final List<String> RequiredOutputAlleles; // limit the alleles with a random distribution written just to save on time and file-size
    public final int Threads;

    private static final String WRITE_RAND_DIST = "write_rand_dist";
    private static final String RANDOM_PEPTIDES_FILE = "random_peptides_file";

    public RandomPeptideConfig(final ConfigBuilder configBuilder)
    {
        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID, "");
        RandomPeptidesFile = configBuilder.getValue(RANDOM_PEPTIDES_FILE);

        String scoreFileDir = checkAddDirSeparator(configBuilder.getValue(SCORE_FILE_DIR));
        String scoreFileId = configBuilder.getValue(SCORE_FILE_ID);

        ScoreDistributionFile = getScoringFilename(configBuilder, scoreFileDir, scoreFileId, FILE_ID_RAND_DIST);
        LikelihoodDistributionFile = getScoringFilename(configBuilder, scoreFileDir, scoreFileId, FILE_ID_LIKELIHOOD_DIST);
        ExpressionLikelihoodDistributionFile = getScoringFilename(configBuilder, scoreFileDir, scoreFileId, FILE_ID_EXP_LIKELIHOOD_DIST);

        WriteRandomDistribution = configBuilder.hasFlag(WRITE_RAND_DIST);

        RequiredOutputAlleles = loadRequiredOutputAlleles(configBuilder.getValue(REQUIRED_OUTPUT_ALLELES));

        Threads = parseThreads(configBuilder);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPathItem(RANDOM_PEPTIDES_FILE, false, "Random peptide file");

        configBuilder.addPathItem(scoreFileConfig(FILE_ID_RAND_DIST), false, "Random peptide distribution file");
        configBuilder.addPathItem(scoreFileConfig(FILE_ID_LIKELIHOOD_DIST), false, "Random likelihood distribution file");
        configBuilder.addPathItem(scoreFileConfig(FILE_ID_EXP_LIKELIHOOD_DIST), false, "Random expression likelihood distribution file");

        configBuilder.addConfigItem(REQUIRED_OUTPUT_ALLELES, "Restricted set of alleles to write to file");
        configBuilder.addFlagItem(WRITE_RAND_DIST, "Write random peptide score distribution");
    }
}
package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.perf.TaskExecutor.parseThreads;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.SCORE_FILE_DIR;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.SCORE_FILE_ID;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.getScoringFilename;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.scoreFileConfig;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_EXP_LIKELIHOOD_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_LIKELIHOOD_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_RAND_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.REQUIRED_OUTPUT_ALLELES;
import static com.hartwig.hmftools.neo.bind.TrainConfig.loadRequiredOutputAlleles;

import java.util.Collections;
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

    public RandomPeptideConfig(
            final String randomPeptidesFile, final String scoreDistributionFile, final String likelihoodDistributionFile,
            final String expressionLikelihoodDistributionFile, final String outputDir, final String outputId,
            final boolean writeRandomDistribution, final List<String> requiredOutputAlleles, final int threads)
    {
        RandomPeptidesFile = randomPeptidesFile;
        ScoreDistributionFile = scoreDistributionFile;
        LikelihoodDistributionFile = likelihoodDistributionFile;
        ExpressionLikelihoodDistributionFile = expressionLikelihoodDistributionFile;
        OutputDir = outputDir;
        OutputId = outputId;
        WriteRandomDistribution = writeRandomDistribution;
        RequiredOutputAlleles = requiredOutputAlleles;
        Threads = threads;
    }

    public static RandomPeptideConfig forReading(final ConfigBuilder configBuilder)
    {
        String scoreFileDir = checkAddDirSeparator(configBuilder.getValue(SCORE_FILE_DIR));
        String scoreFileId = configBuilder.getValue(SCORE_FILE_ID);

        String scoreDistributionFile = getScoringFilename(configBuilder, scoreFileDir, scoreFileId, FILE_ID_RAND_DIST);
        String likelihoodDistributionFile = getScoringFilename(configBuilder, scoreFileDir, scoreFileId, FILE_ID_LIKELIHOOD_DIST);
        String expressionLikelihoodDistributionFile = getScoringFilename(configBuilder, scoreFileDir, scoreFileId, FILE_ID_EXP_LIKELIHOOD_DIST);

        return new RandomPeptideConfig(
                null, scoreDistributionFile, likelihoodDistributionFile, expressionLikelihoodDistributionFile,
                null, null, false, Collections.emptyList(), 0);
    }

    public static void addConfigForReading(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(scoreFileConfig(FILE_ID_RAND_DIST), false, "Random peptide distribution file");
        configBuilder.addPath(scoreFileConfig(FILE_ID_LIKELIHOOD_DIST), false, "Random likelihood distribution file");
        configBuilder.addPath(scoreFileConfig(FILE_ID_EXP_LIKELIHOOD_DIST), false, "Random expression likelihood distribution file");
    }

    public static RandomPeptideConfig forWriting(final ConfigBuilder configBuilder)
    {
        String outputDir = parseOutputDir(configBuilder);
        String outputId = configBuilder.getValue(OUTPUT_ID, "");
        String randomPeptidesFile = configBuilder.getValue(RANDOM_PEPTIDES_FILE);

        boolean writeRandomDistribution = configBuilder.hasFlag(WRITE_RAND_DIST);

        List<String> requiredOutputAlleles = loadRequiredOutputAlleles(configBuilder.getValue(REQUIRED_OUTPUT_ALLELES));
        int threads = parseThreads(configBuilder);

        return new RandomPeptideConfig(
                randomPeptidesFile, null, null, null,
                outputDir, outputId, writeRandomDistribution, requiredOutputAlleles, threads);
    }

    public static void addConfigForWriting(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(RANDOM_PEPTIDES_FILE, false, "Random peptide file");
        configBuilder.addConfigItem(REQUIRED_OUTPUT_ALLELES, "Restricted set of alleles to write to file");
        configBuilder.addFlag(WRITE_RAND_DIST, "Write random peptide score distribution");
    }
}
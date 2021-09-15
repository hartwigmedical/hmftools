package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.OUTPUT_ID;
import static com.hartwig.hmftools.neo.NeoCommon.THREADS;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.SCORE_FILE_DIR;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.SCORE_FILE_ID;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.getScoringFilename;
import static com.hartwig.hmftools.neo.bind.ScoreConfig.scoreFileConfig;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_EXP_LIKELIHOOD_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_LIKELIHOOD_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.FILE_ID_RAND_DIST;
import static com.hartwig.hmftools.neo.bind.TrainConfig.REQUIRED_OUTPUT_ALLELES;
import static com.hartwig.hmftools.neo.bind.TrainConfig.loadRequiredOutputAlleles;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;

import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

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

    public RandomPeptideConfig(final CommandLine cmd)
    {
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID, "");
        RandomPeptidesFile = cmd.getOptionValue(RANDOM_PEPTIDES_FILE);

        String scoreFileDir = cmd.getOptionValue(SCORE_FILE_DIR);
        String scoreFileId = cmd.getOptionValue(SCORE_FILE_ID);

        ScoreDistributionFile = getScoringFilename(cmd, scoreFileDir, scoreFileId, FILE_ID_RAND_DIST);
        LikelihoodDistributionFile = getScoringFilename(cmd, scoreFileDir, scoreFileId, FILE_ID_LIKELIHOOD_DIST);
        ExpressionLikelihoodDistributionFile = getScoringFilename(cmd, scoreFileDir, scoreFileId, FILE_ID_EXP_LIKELIHOOD_DIST);

        WriteRandomDistribution = cmd.hasOption(WRITE_RAND_DIST);

        RequiredOutputAlleles = loadRequiredOutputAlleles(cmd.getOptionValue(REQUIRED_OUTPUT_ALLELES));

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(RANDOM_PEPTIDES_FILE, true, "Random peptide file");

        options.addOption(scoreFileConfig(FILE_ID_RAND_DIST), true, "Random peptide distribution file");
        options.addOption(scoreFileConfig(FILE_ID_LIKELIHOOD_DIST), true, "Random likelihood distribution file");
        options.addOption(scoreFileConfig(FILE_ID_EXP_LIKELIHOOD_DIST), true, "Random expression likelihood distribution file");

        options.addOption(HLA_DEFINITIONS_FILE, true, "HLA allele definitions file");
        options.addOption(REQUIRED_OUTPUT_ALLELES, true, "Restricted set of alleles to write to file");
        options.addOption(WRITE_RAND_DIST, false, "Write random peptide score distribution");
        options.addOption(THREADS, true, "Thread count");
    }
}
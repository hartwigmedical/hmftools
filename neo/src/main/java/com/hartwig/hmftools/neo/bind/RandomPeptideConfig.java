package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.bind.BinderConfig.FILE_ID_LIKELIHOOD_RAND_DIST;
import static com.hartwig.hmftools.neo.bind.BinderConfig.FILE_ID_RAND_DIST;
import static com.hartwig.hmftools.neo.bind.BinderConfig.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.BinderConfig.REQUIRED_OUTPUT_ALLELES;
import static com.hartwig.hmftools.neo.bind.BinderConfig.SCORE_FILE_DIR;
import static com.hartwig.hmftools.neo.bind.BinderConfig.SCORE_FILE_ID;
import static com.hartwig.hmftools.neo.bind.BinderConfig.THREADS;
import static com.hartwig.hmftools.neo.bind.BinderConfig.getScoringFilename;
import static com.hartwig.hmftools.neo.bind.BinderConfig.loadRequiredOutputAlleles;
import static com.hartwig.hmftools.neo.bind.BinderConfig.scoreFileConfig;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class RandomPeptideConfig
{
    public final String RandomPeptidesFile; // list of random peptides from the proteome
    public final String RandomPeptideDistributionFile; // internally scored and ranked distribuion of random peptides
    public final String RandomLikelihoodDistributionFile;
    public final String OutputDir;
    public final String OutputId;

    public final boolean WriteRandomDistribution;
    public final List<String> RequiredOutputAlleles; // limit the alleles with a random distribution written just to save on time and file-size
    public final List<String> WriteRandomScoreAlleles;
    public final int Threads;

    private static final String WRITE_RAND_DIST = "write_rand_dist";
    private static final String RANDOM_PEPTIDES_FILE = "random_peptides_file";
    private static final String WRITE_RAND_SCORE_ALLELES = "write_rand_scores";

    // private static final String RANDOM_PEPTIDE_DIST_FILE = "rand_dist_file";
    // private static final String RANDOM_LIKELIHOOD_DIST_FILE = "likelihood_rand_dist_file";

    public RandomPeptideConfig(final CommandLine cmd)
    {
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID, "");
        RandomPeptidesFile = cmd.getOptionValue(RANDOM_PEPTIDES_FILE);

        String scoreFileDir = cmd.getOptionValue(SCORE_FILE_DIR);
        String scoreFileId = cmd.getOptionValue(SCORE_FILE_ID);

        RandomPeptideDistributionFile = getScoringFilename(cmd, scoreFileDir, scoreFileId, FILE_ID_RAND_DIST);
        RandomLikelihoodDistributionFile = getScoringFilename(cmd, scoreFileDir, scoreFileId, FILE_ID_LIKELIHOOD_RAND_DIST);

        WriteRandomDistribution = cmd.hasOption(WRITE_RAND_DIST);

        RequiredOutputAlleles = loadRequiredOutputAlleles(cmd.getOptionValue(REQUIRED_OUTPUT_ALLELES));

        WriteRandomScoreAlleles = Lists.newArrayList();

        if(cmd.hasOption(WRITE_RAND_SCORE_ALLELES))
        {
            WriteRandomScoreAlleles.addAll(Arrays.stream(cmd.getOptionValue(WRITE_RAND_SCORE_ALLELES).split(ITEM_DELIM)).collect(Collectors.toList()));
        }

        Threads = Integer.parseInt(cmd.getOptionValue(THREADS, "0"));
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(RANDOM_PEPTIDES_FILE, true, "Random peptide file");

        options.addOption(scoreFileConfig(FILE_ID_RAND_DIST), true, "Random peptide distribution file");
        options.addOption(scoreFileConfig(FILE_ID_LIKELIHOOD_RAND_DIST), true, "Random likelihood distribution file");

        options.addOption(HLA_DEFINITIONS_FILE, true, "HLA allele definitions file");
        options.addOption(REQUIRED_OUTPUT_ALLELES, true, "Restricted set of alleles to write to file");
        options.addOption(WRITE_RAND_DIST, false, "Write random peptide score distribution");
        options.addOption(WRITE_RAND_SCORE_ALLELES, true, "List of alleles to write random score data");
        options.addOption(THREADS, true, "Thread count");
    }
}
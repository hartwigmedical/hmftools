package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.bind.BinderConfig.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.BinderConfig.REQUIRED_OUTPUT_ALLELES;
import static com.hartwig.hmftools.neo.bind.BinderConfig.loadRequiredOutputAlleles;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;

import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class RandomPeptideConfig
{
    public final String RandomPeptidesFile; // list of random peptides from the proteome
    public final String RandomPeptideDistributionFile; // internally scored and ranked distribuion of random peptides
    public final String RandomLikelihoodDistributionFile;
    public final String OutputDir;
    public final String OutputId;

    public final boolean WriteRandomDistribution;
    public final List<String> RequiredOutputAlleles; // limit the alleles with a random distribution written just to save on time and file-size

    private static final String WRITE_RAND_DIST = "write_rand_dist";
    private static final String RANDOM_PEPTIDES_FILE = "random_peptides_file";
    private static final String RANDOM_PEPTIDE_DIST_FILE = "random_peptide_dist_file";
    private static final String RANDOM_LIKELIHOOD_DIST_FILE = "random_likelihood_dist_file";

    public RandomPeptideConfig(final CommandLine cmd)
    {
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID, "");
        RandomPeptidesFile = cmd.getOptionValue(RANDOM_PEPTIDES_FILE);
        RandomPeptideDistributionFile = cmd.getOptionValue(RANDOM_PEPTIDE_DIST_FILE);
        RandomLikelihoodDistributionFile = cmd.getOptionValue(RANDOM_LIKELIHOOD_DIST_FILE);
        WriteRandomDistribution = cmd.hasOption(WRITE_RAND_DIST);

        RequiredOutputAlleles = loadRequiredOutputAlleles(cmd.getOptionValue(REQUIRED_OUTPUT_ALLELES));
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(RANDOM_PEPTIDES_FILE, true, "Random peptide file");
        options.addOption(RANDOM_PEPTIDE_DIST_FILE, true, "Random peptide distribution file");
        options.addOption(RANDOM_LIKELIHOOD_DIST_FILE, true, "Random likelihood distribution file");
        options.addOption(HLA_DEFINITIONS_FILE, true, "HLA allele definitions file");
        options.addOption(REQUIRED_OUTPUT_ALLELES, true, "Restricted set of alleles to write to file");
        options.addOption(WRITE_RAND_DIST, false, "Write random peptide score distribution");
    }
}
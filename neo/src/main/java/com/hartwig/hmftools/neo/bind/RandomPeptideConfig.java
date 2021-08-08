package com.hartwig.hmftools.neo.bind;

import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BinderConfig.OUTPUT_ID;
import static com.hartwig.hmftools.neo.bind.HlaSequences.HLA_DEFINITIONS_FILE;
import static com.hartwig.hmftools.neo.bind.HlaSequences.POSITION_HLA_AA_FILE;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.compress.utils.Lists;

public class RandomPeptideConfig
{
    public final String RandomPeptidesFile; // list of random peptides from the proteome
    public final String RandomPeptideDistributionFile; // internally scored and ranked distribuion of random peptides
    public final String OutputDir;
    public final String OutputId;

    public final boolean WriteRandomDistribution;
    public final List<String> AllelesToWrite;

    private static final String WRITE_RAND_DIST = "write_rand_dist";
    private static final String RANDOM_PEPTIDES_FILE = "random_peptides_file";
    private static final String RANDOM_PEPTIDE_DIST_FILE = "random_peptide_dist_file";
    private static final String ALLELES_TO_WRITE = "random_peptide_alleles";

    public RandomPeptideConfig(final CommandLine cmd)
    {
        OutputDir = parseOutputDir(cmd);
        OutputId = cmd.getOptionValue(OUTPUT_ID, "");
        RandomPeptidesFile = cmd.getOptionValue(RANDOM_PEPTIDES_FILE);
        RandomPeptideDistributionFile = cmd.getOptionValue(RANDOM_PEPTIDE_DIST_FILE);
        WriteRandomDistribution = cmd.hasOption(WRITE_RAND_DIST);

        AllelesToWrite = Lists.newArrayList();

        if(cmd.hasOption(ALLELES_TO_WRITE))
        {
            String requiredAlleles = cmd.getOptionValue(ALLELES_TO_WRITE);

            try
            {
                Files.readAllLines(Paths.get(requiredAlleles)).stream()
                        .filter(x -> !x.equals("Allele")).forEach(x -> AllelesToWrite.add(x));

                NE_LOGGER.info("loaded {} required alleles from file({})", AllelesToWrite.size(), requiredAlleles);
            }
            catch(IOException e)
            {
                NE_LOGGER.error("failed to load random peptide alleles to write from file({}): {}", requiredAlleles, e.toString());
            }
        }
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(RANDOM_PEPTIDES_FILE, true, "Random peptide file");
        options.addOption(RANDOM_PEPTIDE_DIST_FILE, true, "Random peptide distribution file");
        options.addOption(HLA_DEFINITIONS_FILE, true, "HLA allele definitions file");
        options.addOption(POSITION_HLA_AA_FILE, true, "Position HLA allele amino acid mapping file");
        options.addOption(ALLELES_TO_WRITE, true, "Restricted set of alleles to write to file");
        options.addOption(WRITE_RAND_DIST, false, "Write random peptide score distribution");
    }
}
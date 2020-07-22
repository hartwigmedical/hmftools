package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.common.amber.AmberPatient;
import com.hartwig.hmftools.common.amber.AmberPatientFactory;
import com.hartwig.hmftools.common.amber.AmberSample;
import com.hartwig.hmftools.common.amber.AmberSampleFactory;
import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.AmberSiteFactory;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthFactory;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.vcf.VCFFileReader;

public class LoadAmberDataV2 {

    private static final Logger LOGGER = LogManager.getLogger(LoadAmberDataV2.class);

    private static final int DEFAULT_MIN_DEPTH = 10;
    private static final double DEFAULT_MIN_HET_AF_PERCENTAGE = 0.4;
    private static final double DEFAULT_MAX_HET_AF_PERCENTAGE = 0.65;

    private static final String SAMPLE = "sample";
    private static final String AMBER_SNP_VCF = "amber_snp_vcf";
    private static final String AMBER_SNP_FILTER_VCF = "mapping_loci_vcf";
    private static final String BED_FILE = "bed";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);

        String tumorSample = cmd.getOptionValue(SAMPLE);
        String amberSnpPath = cmd.getOptionValue(AMBER_SNP_VCF);
        String mappingLoci = cmd.getOptionValue(AMBER_SNP_FILTER_VCF);

        LOGGER.info("Loading mapping loci: {}", mappingLoci);
        final ListMultimap<Chromosome, AmberSite> mappingSites = AmberSiteFactory.sites(mappingLoci);

        final GenomePositionSelector<AmberSite> selector = GenomePositionSelectorFactory.create(mappingSites);

        try (final DatabaseAccess dbAccess = databaseAccess(cmd);
                final VCFFileReader fileReader = new VCFFileReader(new File(amberSnpPath), false)) {

            LOGGER.info("Loading vcf snp data: {}", amberSnpPath);
            final List<BaseDepth> baseDepths = fileReader.iterator()
                    .stream()
                    .map(BaseDepthFactory::fromVariantContext)
                    .filter(x -> selector.select(x).isPresent())
                    .collect(Collectors.toList());

            final AmberSampleFactory amberSampleFactory =
                    new AmberSampleFactory(DEFAULT_MIN_DEPTH, DEFAULT_MIN_HET_AF_PERCENTAGE, DEFAULT_MAX_HET_AF_PERCENTAGE);
            final AmberSample sample = amberSampleFactory.fromBaseDepth(tumorSample, baseDepths);

            LOGGER.info("Loading existing sample data");
            final List<AmberSample> allSamples = dbAccess.readAmberSamples();

            LOGGER.info("Loading existing patient data");
            final List<AmberPatient> allPatients = dbAccess.readAmberPatients();

            LOGGER.info("Finding matches");
            List<AmberPatient> updatedPatients = AmberPatientFactory.create(0.9, sample, allSamples, allPatients);

            LOGGER.info("Writing amber snp data");
            dbAccess.writeAmberSample(sample);

            LOGGER.info("Writing match data");
            dbAccess.writeAmberPatients(tumorSample, updatedPatients);
        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample");
        options.addOption(AMBER_SNP_VCF, true, "Path to the amber directory");
        options.addOption(AMBER_SNP_FILTER_VCF, true, "Path to the amber directory");
        options.addOption(DB_USER, true, "Database user name");
        options.addOption(DB_PASS, true, "Database password");
        options.addOption(DB_URL, true, "Database url");
        options.addOption(BED_FILE, true, "Location of bed file");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String[] args, @NotNull Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull CommandLine cmd) throws SQLException {
        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}

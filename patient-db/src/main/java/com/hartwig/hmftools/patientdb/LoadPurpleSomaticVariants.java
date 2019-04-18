package com.hartwig.hmftools.patientdb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.RefGenomeEnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class LoadPurpleSomaticVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadPurpleSomaticVariants.class);

    private static final String SAMPLE = "sample";
    private static final String VCF = "somatic_vcf";
    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";
    private static final String REF_GENOME = "ref_genome";

    public static void main(@NotNull final String[] args) throws ParseException, IOException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final DatabaseAccess dbAccess = databaseAccess(cmd);

        final String tumorSample = cmd.getOptionValue(SAMPLE);
        final String vcfPath = cmd.getOptionValue(VCF);
        final IndexedFastaSequenceFile sequenceFile = new IndexedFastaSequenceFile(new File(cmd.getOptionValue(REF_GENOME)));

        LOGGER.info("Reading data from {}", vcfPath);
        final SomaticVariantFactory factory = SomaticVariantFactory.unfilteredInstance();
        final List<SomaticVariant> variants = factory.fromVCFFile(tumorSample, vcfPath);

        LOGGER.info("Enriching variants");
        final List<EnrichedSomaticVariant> enrichedVariants = new RefGenomeEnrichedSomaticVariantFactory(sequenceFile).enrich(variants);

        LOGGER.info("Persisting to db");
        dbAccess.writeSomaticVariants(tumorSample, enrichedVariants);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");
        options.addOption(VCF, true, "Path to the structural variant VCF file.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        options.addOption(REF_GENOME, true, "Optional path to the (indexed) ref genome fasta file.");

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static DatabaseAccess databaseAccess(@NotNull final CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }

}

package com.hartwig.hmftools.purple.structural;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.purple.LoadSomaticVariants;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class LoadStructuralVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadSomaticVariants.class);

    private static final String VCF_FILE = "vcf_file";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull final String[] args)
            throws ParseException, IOException, HartwigException, SQLException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);
        final String vcfFileLocation = cmd.getOptionValue(VCF_FILE);
        final DatabaseAccess dbAccess = databaseAccess(cmd);
        final StructuralVariantFactory factory = new StructuralVariantFactory();

        LOGGER.info("Reading VCF File");
        try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                vcfFileLocation, new VCFCodec(), false)) {
            reader.iterator().forEach(factory::addVariantContext);
        }

        LOGGER.info("Persisting variants to database");
        final SampleNames names = new SampleNames(factory.sampleNames());
        dbAccess.writeStructuralVariants(names.tumor(), factory.results());

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(VCF_FILE, true, "Path to the vcf file.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static DatabaseAccess databaseAccess(CommandLine cmd) throws SQLException {
        final String userName = cmd.getOptionValue(DB_USER);
        final String password = cmd.getOptionValue(DB_PASS);
        final String databaseUrl = cmd.getOptionValue(DB_URL);  //e.g. mysql://localhost:port/database";
        final String jdbcUrl = "jdbc:" + databaseUrl;
        return new DatabaseAccess(userName, password, jdbcUrl);
    }
}

package com.hartwig.hmftools.patientdb;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.SomaticVariantStreamWriter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadPurpleSomaticVariants {

    private static final Logger LOGGER = LogManager.getLogger(LoadPurpleSomaticVariants.class);

    private static final String SAMPLE = "sample";
    private static final String REFERENCE = "reference";
    private static final String RNA = "rna";

    private static final String SOMATIC_VCF = "somatic_vcf";

    private static final String DB_USER = "db_user";
    private static final String DB_PASS = "db_pass";
    private static final String DB_URL = "db_url";

    public static void main(@NotNull String[] args) throws ParseException, IOException, SQLException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);
        DatabaseAccess dbAccess = databaseAccess(cmd);

        String tumorSample = cmd.getOptionValue(SAMPLE);
        String referenceSample = cmd.getOptionValue(REFERENCE, null);
        String rnaSample = cmd.getOptionValue(RNA, null);

        String somaticVcf = cmd.getOptionValue(SOMATIC_VCF);

        LOGGER.info("Removing old data of sample {}", tumorSample);
        try (SomaticVariantStreamWriter somaticWriter = dbAccess.somaticVariantWriter(tumorSample)) {
            LOGGER.info("Streaming data from {} to db", somaticVcf);
            new SomaticVariantFactory().fromVCFFile(tumorSample, referenceSample, rnaSample, somaticVcf, somaticWriter);
        }

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the tumor sample. This should correspond to the value used in PURPLE.");
        options.addOption(REFERENCE, true, "Optional name of the reference sample. This should correspond to the value used in PURPLE.");
        options.addOption(RNA, true, "Optional name of the rna sample. This should correspond to the value used in PURPLE.");
        options.addOption(SOMATIC_VCF, true, "Path to the PURPLE somatic variant VCF file.");
        options.addOption(DB_USER, true, "Database user name.");
        options.addOption(DB_PASS, true, "Database password.");
        options.addOption(DB_URL, true, "Database url.");

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

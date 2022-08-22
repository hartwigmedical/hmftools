package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.sql.SQLException;

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class UpdateSnpCheckStatus {

    private static final Logger LOGGER = LogManager.getLogger(UpdateSnpCheckStatus.class);

    private static final String SAMPLE = "sample";
    private static final String ISOLATION_BARCODE = "isolation_barcode";
    private static final String IS_PASS = "is_pass";

    private static final String HAS_PASSED = "true";
    private static final String HAS_FAILED = "false";

    public static void main(@NotNull String[] args) throws ParseException, SQLException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);
        String isolationBarcode = cmd.getOptionValue(ISOLATION_BARCODE);
        String isPassString = cmd.getOptionValue(IS_PASS);

        if (Utils.anyNull(sample, isolationBarcode, isPassString) || (!isPassString.equals(HAS_PASSED)
                && !isPassString.equals(HAS_FAILED))) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Update SnpCheck Status", options);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(cmd);

        assert isPassString.equals(HAS_PASSED) || isPassString.equals(HAS_FAILED);
        boolean isPass = isPassString.equals(HAS_PASSED);
        LOGGER.info("Updating SnpCheck for {} to '{}'", sample, isPass);
        dbWriter.writeSnpCheck(sample, isolationBarcode, isPass);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "Sample for which we are going to update snpcheck status");
        options.addOption(ISOLATION_BARCODE, true, "Isolation barcode for which we are going to update snpcheck status");
        options.addOption(IS_PASS, true, "Either pass 'true' or 'false' to set snpcheck pass status");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}

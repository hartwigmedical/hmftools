package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.patientdb.bqr.BQREntry;
import com.hartwig.hmftools.patientdb.bqr.BQRFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadSageBQRData {

    private static final Logger LOGGER = LogManager.getLogger(LoadSageBQRData.class);

    private static final String SAMPLE = "sample";

    private static final String REF_SAMPLE_BQR_TSV = "ref_sample_bqr_tsv";
    private static final String TUMOR_SAMPLE_BQR_TSV = "tumor_sample_bqr_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);
        String refBQRTsv = cmd.getOptionValue(REF_SAMPLE_BQR_TSV);
        String tumorBQRTsv = cmd.getOptionValue(TUMOR_SAMPLE_BQR_TSV);

        if (Utils.anyNull(sample, refBQRTsv, tumorBQRTsv) || !new File(refBQRTsv).exists() || !new File(tumorBQRTsv).exists()) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load Sage BQR Data", options);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(cmd);

        LOGGER.info("Reading BQR ref data for {} from {}", sample, refBQRTsv);
        List<BQREntry> refEntries = BQRFile.read(refBQRTsv);

        LOGGER.info("Reading BQR tumor data for {} from {}", sample, tumorBQRTsv);
        List<BQREntry> tumorEntries = BQRFile.read(tumorBQRTsv);

        dbWriter.writeBQR(sample, refEntries, tumorEntries);
        LOGGER.info("Done!");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        options.addOption(SAMPLE, true, "The tumor sample.");
        options.addOption(REF_SAMPLE_BQR_TSV, true, "Path towards the ref sample BQR tsv.");
        options.addOption(TUMOR_SAMPLE_BQR_TSV, true, "Path towards the tumor sample BQR tsv.");
        addDatabaseCmdLineArgs(options);
        return options;
    }
}

package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.common.virusbreakend.VirusBreakend;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadVirusBreakendData {

    private static final Logger LOGGER = LogManager.getLogger(LoadVirusBreakendData.class);

    private static final String SAMPLE = "sample";
    private static final String VIRUS_BREAKEND_TSV = "virus_breakend_tsv";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);
        String viusBreakendTsv = cmd.getOptionValue(VIRUS_BREAKEND_TSV);

        if (Utils.anyNull(sample, viusBreakendTsv)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load VIRUSBreakend Data", options);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(cmd);

        LOGGER.info("Reading virus breakend TSV {}", viusBreakendTsv);
        List<VirusBreakend> virusBreakends = VirusBreakendFile.read(viusBreakendTsv);
        LOGGER.info(" Read {} virus breakends", virusBreakends.size());

        LOGGER.info("Writing virus breakends into database for {}", sample);
        dbWriter.writeVirusBreakend(sample, virusBreakends);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(SAMPLE, true, "Sample for which we are going to load the virus breakends");
        options.addOption(VIRUS_BREAKEND_TSV, true, "Path towards the virus breakend TSV file");

        addDatabaseCmdLineArgs(options);

        return options;
    }

}

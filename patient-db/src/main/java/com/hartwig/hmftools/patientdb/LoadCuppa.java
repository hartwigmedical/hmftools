package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.databaseAccess;

import java.io.IOException;
import java.sql.SQLException;

import com.hartwig.hmftools.common.cuppa.MolecularTissueOrginData;
import com.hartwig.hmftools.common.cuppa.MolecularTissueOriginFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class LoadCuppa {

    private static final Logger LOGGER = LogManager.getLogger(LoadCuppa.class);

    private static final String SAMPLE = "sample";
    private static final String CUPPA_CONCLUSION_TXT = "cuppa_conclusion_txt";

    public static void main(@NotNull String[] args) throws ParseException, SQLException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        String sample = cmd.getOptionValue(SAMPLE);
        String cuppaConclusionTxt = cmd.getOptionValue(CUPPA_CONCLUSION_TXT);

        if (Utils.anyNull(sample, cuppaConclusionTxt)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Load CUPPA Data", options);
            System.exit(1);
        }

        DatabaseAccess dbWriter = databaseAccess(cmd);

        LOGGER.info("Reading CUPPA conclusion file {}", cuppaConclusionTxt);
        MolecularTissueOrginData molecularTissueOrigins = MolecularTissueOriginFile.read(cuppaConclusionTxt);
        LOGGER.info(" Read '{}' as CUPPA conclusion result", molecularTissueOrigins.conclusion());

        LOGGER.info("Writing CUPPA into database for {}", sample);
        dbWriter.writeCuppa(sample, molecularTissueOrigins);

        LOGGER.info("Complete");
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(SAMPLE, true, "Sample for which we are going to load the CUPPA results");
        options.addOption(CUPPA_CONCLUSION_TXT, true, "Path towards the CUPPA conclusion txt file");

        addDatabaseCmdLineArgs(options);

        return options;
    }
}
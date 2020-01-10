package com.hartwig.hmftools.iclusion;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import com.hartwig.hmftools.iclusion.api.iclusionApi;

public class IclusionImporterApplication {
    private static final Logger LOGGER = LogManager.getLogger(IclusionImporterApplication.class);
    private static final String ICLUSION_LINK = "iclusion_link";
    private static final String ICLUSION_CLIENT_ID = "iclusion_client_id";
    private static final String ICLUSION_CLIENT_SECRET = "iclusion_client_secret";
    private static final String ICLUSION_USERNAME = "iclusion_username";
    private static final String ICLUSION_PASSWORD = "iclusion_password";

    private static final String ICLUSION_OUTPUT_STUDIES_RAW = "iclusion_output_studies_raw";
    private static final String ICLUSION_OUTPUT_STUDIES_PROCESSED = "iclusion_output_studies_processed";

    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (!validInputForBaseReport(cmd)) {
            printUsageAndExit(options);
        }

        String iClusionLink = cmd.getOptionValue(ICLUSION_LINK);
        String iClusionClientId = cmd.getOptionValue(ICLUSION_CLIENT_ID);
        String iClusionClientSecret = cmd.getOptionValue(ICLUSION_CLIENT_SECRET);
        String iClusionUsername = cmd.getOptionValue(ICLUSION_USERNAME);
        String iClusionPassword = cmd.getOptionValue(ICLUSION_PASSWORD);
        String iClusionOutputStudiesRaw = cmd.getOptionValue(ICLUSION_OUTPUT_STUDIES_RAW);
        String iClusionOutputStudiesProcessed = cmd.getOptionValue(ICLUSION_OUTPUT_STUDIES_PROCESSED);

        String token = iclusionApi.connectWithIclusionApi(iClusionLink,
                iClusionClientId,
                iClusionClientSecret,
                iClusionUsername,
                iClusionPassword);

        LOGGER.info("Reading iclusion study details.....");
        LOGGER.info("Queried and filtered {} studies from iclusion API", "size study");

        writeIclusionOutputStudiesRawToTSVFile(iClusionOutputStudiesRaw);
        writeIclusionOutputStudiesProcessedToTSVFile(iClusionOutputStudiesProcessed);

        LOGGER.info("Iclusion importer is finished!");

    }

    private static void writeIclusionOutputStudiesRawToTSVFile(@NotNull String iClusionOutputStudiesRaw) throws IOException {
        LOGGER.info("Writing iClusion output raw to file {}", iClusionOutputStudiesRaw);
        BufferedWriter writer = new BufferedWriter(new FileWriter(iClusionOutputStudiesRaw, true));
        writer.write("gene"); //TODO write real data from iclusion
        writer.write(""); //TODO write real data from iclusion
        writer.close();
    }

    private static void writeIclusionOutputStudiesProcessedToTSVFile(@NotNull String iClusionOutputStudiesProcessed) throws IOException {
        LOGGER.info("Writing iClusion output processed to file {}", iClusionOutputStudiesProcessed);
        BufferedWriter writer = new BufferedWriter(new FileWriter(iClusionOutputStudiesProcessed, true));
        writer.write("gene"); //TODO write real data from iclusion
        writer.write(""); //TODO write real data from iclusion
        writer.close();
    }

    private static boolean validInputForBaseReport(@NotNull CommandLine cmd) {
        return valueExists(cmd, ICLUSION_LINK) && valueExists(cmd, ICLUSION_CLIENT_ID) && valueExists(cmd, ICLUSION_CLIENT_SECRET)
                && valueExists(cmd, ICLUSION_USERNAME) && valueExists(cmd, ICLUSION_PASSWORD);
    }

    private static boolean valueExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);
        if (value == null) {
            LOGGER.warn(param + " has to be provided");
            return false;
        }
        return true;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(ICLUSION_LINK, true, "iClusion link");
        options.addOption(ICLUSION_CLIENT_ID, true, "iClusion client id");
        options.addOption(ICLUSION_CLIENT_SECRET, true, "iClusion client secret");
        options.addOption(ICLUSION_USERNAME, true, "iClusion username");
        options.addOption(ICLUSION_PASSWORD, true, "iClusion passsword");
        options.addOption(ICLUSION_OUTPUT_STUDIES_RAW, true, "iClusion output studies raw");
        options.addOption(ICLUSION_OUTPUT_STUDIES_PROCESSED, true, "iClusion output studies processed");

        return options;
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Iclusion-Importer", options);
        System.exit(1);
    }
}

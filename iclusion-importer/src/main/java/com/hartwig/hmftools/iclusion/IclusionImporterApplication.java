package com.hartwig.hmftools.iclusion;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionImporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(IclusionImporterApplication.class);

    private static final String ICLUSION_ENDPOINT = "iclusion_endpoint";
    private static final String ICLUSION_CLIENT_ID = "iclusion_client_id";
    private static final String ICLUSION_CLIENT_SECRET = "iclusion_client_secret";
    private static final String ICLUSION_USERNAME = "iclusion_username";
    private static final String ICLUSION_PASSWORD = "iclusion_password";

    private static final String ICLUSION_OUTPUT_STUDIES_RAW = "iclusion_output_studies_raw";
    private static final String ICLUSION_OUTPUT_STUDIES_PROCESSED = "iclusion_output_studies_processed";

    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        Options options = createBasicOptions();
        CommandLine cmd = createCommandLine(args, options);

        if (!validInputForIclusionConnection(cmd)) {
            printUsageAndExit(options);
        }

        String iClusionEndpoint = cmd.getOptionValue(ICLUSION_ENDPOINT);
        String iClusionClientId = cmd.getOptionValue(ICLUSION_CLIENT_ID);
        String iClusionClientSecret = cmd.getOptionValue(ICLUSION_CLIENT_SECRET);
        String iClusionUsername = cmd.getOptionValue(ICLUSION_USERNAME);
        String iClusionPassword = cmd.getOptionValue(ICLUSION_PASSWORD);

//        String token = IclusionApiWrapper.connectWithIclusionApi(iClusionEndpoint,
//                iClusionClientId,
//                iClusionClientSecret,
//                iClusionUsername,
//                iClusionPassword);

        String iClusionOutputStudiesRaw = cmd.getOptionValue(ICLUSION_OUTPUT_STUDIES_RAW);
        String iClusionOutputStudiesProcessed = cmd.getOptionValue(ICLUSION_OUTPUT_STUDIES_PROCESSED);

        LOGGER.info("Reading iClusion study details.....");
        LOGGER.info("Queried and filtered {} studies from iClusion API", "size study");

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

    private static boolean validInputForIclusionConnection(@NotNull CommandLine cmd) {
        return valueExists(cmd, ICLUSION_ENDPOINT) && valueExists(cmd, ICLUSION_CLIENT_ID) && valueExists(cmd, ICLUSION_CLIENT_SECRET)
                && valueExists(cmd, ICLUSION_USERNAME) && valueExists(cmd, ICLUSION_PASSWORD);
    }

    private static boolean valueExists(@NotNull CommandLine cmd, @NotNull String param) {
        if (!cmd.hasOption(param)) {
            LOGGER.warn("{} has to be provided", param);
            return false;
        }
        return true;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        return new DefaultParser().parse(options, args);
    }

    @NotNull
    private static Options createBasicOptions() {
        Options options = new Options();

        options.addOption(ICLUSION_ENDPOINT, true, "iClusion endpoint URL");
        options.addOption(ICLUSION_CLIENT_ID, true, "iClusion client id");
        options.addOption(ICLUSION_CLIENT_SECRET, true, "iClusion client secret");
        options.addOption(ICLUSION_USERNAME, true, "iClusion username");
        options.addOption(ICLUSION_PASSWORD, true, "iClusion password");
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

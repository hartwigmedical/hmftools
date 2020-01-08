package com.hartwig.hmftools.iclusion;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URL;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import sun.net.www.http.HttpClient;

public class IclusionImporterApplication {
   // private static final Logger LOGGER = LogManager.getLogger(IclusionImporterApplication.class);
    private static final String ICLUSION_LINK = "iclusion_link";
    private static final String ICLUSION_CLIENT_ID = "iclusion_client_id";
    private static final String ICLUSION_CLIENT_SECRET = "iclusion_client_secret";
    private static final String ICLUSION_USERNAME = "iclusion_username";
    private static final String ICLUSION_PASSWORD = "iclusion_password";

    private static final String ICLUSION_OUTPUT_STUDIES = "iclusion_output_studies";

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
        String iClusionOutputStudies = cmd.getOptionValue(ICLUSION_OUTPUT_STUDIES);

    //    LOGGER.info("Connecting with iclusion API on {}", iClusionLink);

        connectWithIclusionApi(iClusionLink,
                iClusionClientId,
                iClusionClientSecret,
                iClusionUsername,
                iClusionPassword,
                iClusionOutputStudies);

        //  LOGGER.info("Reading iclusion study details.....");
        //  LOGGER.info("Queried and filtered {} studies from iclusion API", "size study");

        writeIclusionOutputStudiesToTSVFile(iClusionOutputStudies);
        //  LOGGER.info("Iclusion importer is finished!");

    }

    private static void connectWithIclusionApi(@NotNull String iClusionLink, @NotNull String iClusionClientId,
            @NotNull String iClusionClientSecret, @NotNull String iClusionUsername, @NotNull String iClusionPassword,
            @NotNull String iClusionOutputStudies) throws IOException {

        URL url = new URL(iClusionLink); // url iclusion
        BufferedWriter writer = new BufferedWriter(new FileWriter(iClusionOutputStudies));
        HttpURLConnection connection = (HttpURLConnection) url.openConnection();

//        connection.setReadTimeout(1000);
//        connection.setConnectTimeout(1000);
//
//        connection.setDoInput(true);
//        connection.setDoOutput(true);
//        connection.setRequestMethod("POST");
//        connection.setRequestProperty("Accept", "application/json");
//        connection.setRequestProperty("content-type", "multipart/form-data");
//        connection.setRequestProperty("grant_type", "password");
//        connection.setRequestProperty("client_id", iClusionClientId);
//        connection.setRequestProperty("client_secret", iClusionClientSecret);
//        connection.setRequestProperty("username", iClusionUsername);
//        connection.setRequestProperty("password", iClusionPassword);
//        //  connection.setRequestProperty("Authorization", "Bearer " + "");
//        connection.connect();

        writer.write(connection.toString() + "\n");
        writer.write(connection.getRequestMethod() + "\n");
        writer.write(connection.getResponseMessage() + "\n");
        writer.write(connection.getResponseCode() + "\n");
        writer.write(connection.getHeaderFields() + "\n");

        writer.close();

    }

    private static void writeIclusionOutputStudiesToTSVFile(@NotNull String iClusionOutputStudies) throws IOException {
        //  LOGGER.info("Writing iClusion output to file {}", iClusionOutputStudies);
        BufferedWriter writer = new BufferedWriter(new FileWriter(iClusionOutputStudies, true));
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
            //   LOGGER.warn(param + " has to be provided");
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
        options.addOption(ICLUSION_OUTPUT_STUDIES, true, "iClusion output studies");
        return options;
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Iclusion-Importer", options);
        System.exit(1);
    }
}

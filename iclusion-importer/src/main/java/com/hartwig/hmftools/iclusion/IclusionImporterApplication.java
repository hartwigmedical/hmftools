package com.hartwig.hmftools.iclusion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class IclusionImporterApplication {
    private static final Logger LOGGER = LogManager.getLogger(IclusionImporterApplication.class);
    private static final String ICLUSION_LINK = "iclusion_link";
    private static final String ICLUSION_CLIENT_ID = "iclusion_client_id";
    private static final String ICLUSION_CLIENT_SECRET = "iclusion_client_secret";
    private static final String ICLUSION_USERNAME = "iclusion_username";
    private static final String ICLUSION_PASSWORD = "iclusion_password";

    public static void main(@NotNull final String[] args) throws ParseException {
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

        LOGGER.info("Connecting with iclusion API on {}", iClusionLink);

        connectWithIclusionApi(iClusionLink, iClusionClientId, iClusionClientSecret, iClusionUsername, iClusionPassword);

        LOGGER.info("Reading iclusion study details.....");
        LOGGER.info("Queried and filtered {} studies from iclusion API", "size study");


        //        val iclusionApi = IclusionApiWrapper(iclusionEndpoint, cmd.getOptionValue(ICLUSION_CLIENT_ID),
        //                cmd.getOptionValue(ICLUSION_CLIENT_SECRET), cmd.getOptionValue(ICLUSION_USER), cmd.getOptionValue(ICLUSION_PASSWORD))
        //        logger.info("Reading iclusion study details...")
        //        val iclusionStudies = iclusionApi.studyDetails()
        //        iclusionApi.close()
        //
        //        iclusionStudies.forEach {
        //            logger.info("iclusion study: $it")
        //        }
        //        logger.info("Queried and filtered ${iclusionStudies.size} studies from iclusion API")
    }

    private static void connectWithIclusionApi(@NotNull String iClusionLink, @NotNull String iClusionClientId,
            @NotNull String iClusionClientSecret, @NotNull String iClusionUsername, @NotNull String iClusionPassword) {


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
        return options;
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Iclusion-Importer", options);
        System.exit(1);
    }
}

package com.hartwig.hmftools.iclusion;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.iclusion.api.IclusionApiMain;
import com.hartwig.hmftools.iclusion.data.IclusionTrial;
import com.hartwig.hmftools.iclusion.io.IclusionTrialFile;
import com.hartwig.hmftools.iclusion.qc.IclusionTrialChecker;

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

    private static final String ICLUSION_TRIAL_TSV = "iclusion_trial_tsv";

    private static final String VERSION = IclusionImporterApplication.class.getPackage().getImplementationVersion();

    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        LOGGER.info("Running iClusion Importer v{}", VERSION);

        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        if (!validInputForIclusionImporting(cmd)) {
            printUsageAndExit(options);
        }

        List<IclusionTrial> trials = IclusionApiMain.readIclusionTrials(cmd.getOptionValue(ICLUSION_ENDPOINT), buildCredentials(cmd));

        IclusionTrialChecker.check(trials);

        String iClusionTrialTsv = cmd.getOptionValue(ICLUSION_TRIAL_TSV);
        IclusionTrialFile.write(iClusionTrialTsv, trials);

        LOGGER.info("iClusion Importer has written {} trials to {}", trials.size(), iClusionTrialTsv);
    }

    @NotNull
    private static IclusionCredentials buildCredentials(@NotNull CommandLine cmd) {
        return ImmutableIclusionCredentials.builder()
                .clientId(cmd.getOptionValue(ICLUSION_CLIENT_ID))
                .clientSecret(cmd.getOptionValue(ICLUSION_CLIENT_SECRET))
                .username(cmd.getOptionValue(ICLUSION_USERNAME))
                .password(cmd.getOptionValue(ICLUSION_PASSWORD))
                .build();
    }

    private static boolean validInputForIclusionImporting(@NotNull CommandLine cmd) {
        return valueExists(cmd, ICLUSION_ENDPOINT) && valueExists(cmd, ICLUSION_CLIENT_ID) && valueExists(cmd, ICLUSION_CLIENT_SECRET)
                && valueExists(cmd, ICLUSION_USERNAME) && valueExists(cmd, ICLUSION_PASSWORD) && valueExists(cmd, ICLUSION_TRIAL_TSV);
    }

    private static boolean valueExists(@NotNull CommandLine cmd, @NotNull String param) {
        if (!cmd.hasOption(param)) {
            LOGGER.warn("{} has to be provided", param);
            return false;
        }
        return true;
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(ICLUSION_ENDPOINT, true, "iClusion endpoint URL");
        options.addOption(ICLUSION_CLIENT_ID, true, "iClusion client id");
        options.addOption(ICLUSION_CLIENT_SECRET, true, "iClusion client secret");
        options.addOption(ICLUSION_USERNAME, true, "iClusion username");
        options.addOption(ICLUSION_PASSWORD, true, "iClusion password");

        options.addOption(ICLUSION_TRIAL_TSV, true, "iClusion output trial tsv");

        return options;
    }

    private static void printUsageAndExit(@NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Iclusion-Importer", options);
        System.exit(1);
    }
}

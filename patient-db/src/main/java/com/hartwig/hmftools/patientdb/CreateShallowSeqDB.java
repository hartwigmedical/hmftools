package com.hartwig.hmftools.patientdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;

import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.common.purple.qc.PurpleQCStatus;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class CreateShallowSeqDB {

    private static final Logger LOGGER = LogManager.getLogger(CreateShallowSeqDB.class);

    private static final String SAMPLE = "sample";

    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_QC_FILE = "purple_qc_file";

    private static final String SHALLOW_SEQ_CSV = "shallow_seq_csv";


    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (validInputForAnalysedSample(cmd)) {
            PurityContext purityContext = FittedPurityFile.read(cmd.getOptionValue(PURPLE_PURITY_TSV));
            PurpleQC purpleQC = PurpleQCFile.read(cmd.getOptionValue(PURPLE_QC_FILE));

            PurpleQCStatus QCstatus = purpleQC.status();
            FittedPurityStatus status = purityContext.status();
            String purity = Strings.EMPTY;
            if (status == FittedPurityStatus.NO_TUMOR) {
                purity = "below detection threshold";
            } else {
                purity = String.valueOf(purityContext.bestFit().purity());
            }
            String outputStringForFile = cmd.getOptionValue(PURPLE_QC_FILE) + "," + purity;

            //TODO: check if sample already exist in file and remove from file
            appendToTsv(cmd.getOptionValue(SHALLOW_SEQ_CSV), outputStringForFile);

        } else {
            printUsageAndExit(options);
        }

        LOGGER.info("Complete");
    }

    private static void appendToTsv(@NotNull String shallowSeqCsv, @NotNull String stringToAppend) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(shallowSeqCsv, true));
        writer.write(stringToAppend);
        writer.close();
    }

    private static boolean validInputForAnalysedSample(@NotNull CommandLine cmd) {
        return fileExists(cmd, PURPLE_PURITY_TSV) && fileExists(cmd, PURPLE_QC_FILE);
    }

    private static boolean fileExists(@NotNull CommandLine cmd, @NotNull String param) {
        String value = cmd.getOptionValue(param);

        if (value == null || !pathExists(value)) {
            LOGGER.warn(param + " has to be an existing file: " + value);
            return false;
        }

        return true;
    }

    private static boolean pathExists(@NotNull String path) {
        return Files.exists(new File(path).toPath());
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Tumor sample.");

        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_QC_FILE, true, "Path towards the purple qc file.");

        options.addOption(SHALLOW_SEQ_CSV, true, "Path towards output file for the shallow seq db CSV.");


        return options;
    }

    private static void printUsageAndExit(@NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Create Shallow Seq DB", options);
        System.exit(1);
    }

}

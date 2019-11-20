package com.hartwig.hmftools.patientdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityStatus;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.common.purple.qc.PurpleQCStatus;
import com.hartwig.hmftools.patientdb.context.RunContext;
import com.hartwig.hmftools.patientdb.readers.RunsFolderReader;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CreateShallowSeqDB {

    private static final Logger LOGGER = LogManager.getLogger(CreateShallowSeqDB.class);

    private static final String RUNS_DIRECTORY = "runs_dir";

    private static final String PURPLE_PURITY_TSV = "purple_purity_tsv";
    private static final String PURPLE_QC_FILE = "purple_qc_file";

    private static final String SHALLOW_SEQ_CSV = "shallow_seq_csv";

    private static final String PURPLE_DIR = "/purple/";

    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        LOGGER.info("Loading shallow seq runs from {}", cmd.getOptionValue(RUNS_DIRECTORY));
        final List<RunContext> runContexts = loadRunContexts(cmd.getOptionValue(RUNS_DIRECTORY));

        extraxtPurpleFromSample(runContexts,
                cmd.getOptionValue(PURPLE_PURITY_TSV),
                cmd.getOptionValue(PURPLE_QC_FILE),
                cmd.getOptionValue(SHALLOW_SEQ_CSV));

        LOGGER.info("Complete");
    }

    @NotNull
    private static void extraxtPurpleFromSample(@NotNull List<RunContext> runContexts, @NotNull String purplePurityTsv,
            @NotNull String purpleQCTsv, @NotNull String shallowSeqOutputCsv) throws IOException {
        for (RunContext runInfo : runContexts) {
            String tumorSample = runInfo.tumorSample();
            String setName = runInfo.setName();
            String sampleBarcode = "";

            PurityContext purityContext = FittedPurityFile.read(setName + PURPLE_DIR + tumorSample + purplePurityTsv);
            PurpleQC purpleQC = PurpleQCFile.read(setName + PURPLE_DIR + tumorSample + purpleQCTsv);

            boolean QCstatus = purpleQC.status() == PurpleQCStatus.PASS;
            boolean status = purityContext.status() != FittedPurityStatus.NO_TUMOR;
            double purity = purityContext.bestFit().purity();

            String outputStringForFile = sampleBarcode + "," + tumorSample + "," + purity + "," + QCstatus + "," + status;
            appendToTsv(shallowSeqOutputCsv, outputStringForFile);

        }
    }

    @NotNull
    private static List<RunContext> loadRunContexts(@NotNull String runsDirectory) throws IOException {
        final List<RunContext> runContexts = RunsFolderReader.extractRunContexts(new File(runsDirectory));
        LOGGER.info(" Loaded run contexts from {} ({} sets)", runsDirectory, runContexts.size());

        return runContexts;
    }

    private static void appendToTsv(@NotNull String shallowSeqCsv, @NotNull String stringToAppend) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(shallowSeqCsv, true));
        writer.write(stringToAppend);
        writer.close();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    @NotNull
    private static Options createBasicOptions() {
        final Options options = new Options();

        options.addOption(RUNS_DIRECTORY, true, "Path towards the folder containing all shallow seq runs .");

        options.addOption(PURPLE_PURITY_TSV, true, "Path towards the purple purity TSV.");
        options.addOption(PURPLE_QC_FILE, true, "Path towards the purple qc file.");

        options.addOption(SHALLOW_SEQ_CSV, true, "Path towards output file for the shallow seq db CSV.");

        return options;
    }

}

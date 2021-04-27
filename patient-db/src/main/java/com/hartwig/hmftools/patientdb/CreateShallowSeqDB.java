package com.hartwig.hmftools.patientdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.ImmutableLimsShallowSeqData;
import com.hartwig.hmftools.common.lims.LimsShallowSeqData;
import com.hartwig.hmftools.common.purple.CheckPurpleQuality;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.runcontext.RunContext;
import com.hartwig.hmftools.patientdb.clinical.readers.RunsFolderReader;

import org.apache.commons.cli.CommandLine;
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

    private static final String RUNS_DIRECTORY = "runs_dir";

    private static final String SHALLOW_SEQ_TSV = "shallow_seq_tsv";
    private static final String PURPLE_PURITY_P4_TSV = "purple_purity_p4_tsv";
    private static final String PURPLE_PURITY_P5_TSV = "purple_purity_p5_tsv";
    private static final String PURPLE_QC_FILE = "purple_qc_file";
    private static final String PIPELINE_VERSION_FILE = "pipeline_version_file";

    private static final String PURPLE_DIR = "purple";
    private static final String DELIMITER = "\t";

    public static void main(@NotNull String[] args) throws ParseException, IOException {
        Options options = createOptions();
        CommandLine cmd = new DefaultParser().parse(options, args);

        if (!checkInputs(cmd)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-DB - Create Shallow Seq DB ", options);
            System.exit(1);
        }

        LOGGER.info("Loading shallow seq runs from {}", cmd.getOptionValue(RUNS_DIRECTORY));
        List<RunContext> runContexts = loadRunContexts(cmd.getOptionValue(RUNS_DIRECTORY), cmd.getOptionValue(PIPELINE_VERSION_FILE));

        List<LimsShallowSeqData> newShallowSeqEntries = extractNewEntriesForShallowDbFromRunContexts(runContexts,
                cmd.getOptionValue(SHALLOW_SEQ_TSV),
                cmd.getOptionValue(RUNS_DIRECTORY),
                cmd.getOptionValue(PURPLE_QC_FILE),
                cmd.getOptionValue(PURPLE_PURITY_P4_TSV),
                cmd.getOptionValue(PURPLE_PURITY_P5_TSV),
                cmd.getOptionValue(PIPELINE_VERSION_FILE));

        appendToCsv(cmd.getOptionValue(SHALLOW_SEQ_TSV), newShallowSeqEntries);

        LOGGER.info("Complete");
    }

    @NotNull
    private static List<LimsShallowSeqData> extractNewEntriesForShallowDbFromRunContexts(@NotNull List<RunContext> runContexts,
            @NotNull String shallowSeqTsv, @NotNull String runDirsPath, @NotNull String purpleQCFile, @NotNull String purplePurityExtP4,
            @NotNull String purplePurityExtP5, @NotNull String pipelineVersionExt) throws IOException {
        List<LimsShallowSeqData> currentShallowSeqData = read(shallowSeqTsv);

        List<LimsShallowSeqData> shallowSeqDataToAppend = Lists.newArrayList();

        for (RunContext runInfo : runContexts) {
            String tumorSample = runInfo.tumorSample();
            String setPath = runDirsPath + File.separator + runInfo.setName();
            String sampleBarcode = runInfo.tumorBarcodeSample();

            String purplePurityTsvExt;
            if (new File(setPath + File.separator + pipelineVersionExt).exists()) {
                purplePurityTsvExt = purplePurityExtP5;
            } else {
                purplePurityTsvExt = purplePurityExtP4;
            }
            String fullPurplePurityTsvPath = setPath + File.separator + PURPLE_DIR + File.separator + tumorSample + purplePurityTsvExt;
            String fullPurpleQCFilePath = setPath + File.separator + PURPLE_DIR + File.separator + tumorSample + purpleQCFile;

            PurityContext purityContext = PurityContextFile.readWithQC(fullPurpleQCFilePath, fullPurplePurityTsvPath);
            PurpleQC purpleQC = purityContext.qc();

            boolean hasReliableQuality = purpleQC.pass();
            boolean hasReliablePurity = CheckPurpleQuality.checkHasReliablePurity(purityContext);
            String purity = new DecimalFormat("0.00").format(purityContext.bestFit().purity());

            boolean inFile = false;
            for (LimsShallowSeqData sample : currentShallowSeqData) {
                if (sample.sampleBarcode().equals(sampleBarcode)) {
                    LOGGER.warn("Sample barcode is already present in file. Skipping set: {} with sample barcode: {} for"
                            + " writing to shallow seq db!", setPath, sampleBarcode);
                    inFile = true;
                }
            }
            if (!inFile && !sampleBarcode.equals(Strings.EMPTY)) {
                shallowSeqDataToAppend.add(ImmutableLimsShallowSeqData.builder()
                        .sampleBarcode(sampleBarcode)
                        .sampleId(tumorSample)
                        .purityShallowSeq(purity)
                        .hasReliableQuality(hasReliableQuality)
                        .hasReliablePurity(hasReliablePurity)
                        .build());
                LOGGER.info("Set: {} is added to shallow list!", setPath);
            }
        }
        return shallowSeqDataToAppend;
    }

    @NotNull
    private static List<LimsShallowSeqData> read(@NotNull String shallowSeqTsv) throws IOException {
        List<String> linesShallowDB = Files.readAllLines(new File(shallowSeqTsv).toPath());
        List<LimsShallowSeqData> shallowSeqDataList = Lists.newArrayList();
        for (String line : linesShallowDB.subList(1, linesShallowDB.size())) {
            String[] values = line.split(DELIMITER);

            shallowSeqDataList.add(ImmutableLimsShallowSeqData.builder()
                    .sampleBarcode(values[0])
                    .sampleId(values[1])
                    .purityShallowSeq(values[2])
                    .hasReliableQuality(Boolean.parseBoolean(values[3]))
                    .hasReliablePurity(Boolean.parseBoolean(values[4]))
                    .build());
        }
        return shallowSeqDataList;
    }

    @NotNull
    private static List<RunContext> loadRunContexts(@NotNull String runsDirectory, @NotNull String pipelineVersionFile) throws IOException {
        List<RunContext> runContexts = RunsFolderReader.extractRunContexts(new File(runsDirectory), pipelineVersionFile);
        LOGGER.info(" Loaded run contexts from {} ({} sets)", runsDirectory, runContexts.size());

        return runContexts;
    }

    private static void appendToCsv(@NotNull String shallowSeqCsv, @NotNull List<LimsShallowSeqData> shallowSeqDataToAppend)
            throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(shallowSeqCsv, true));
        for (LimsShallowSeqData dataToAppend : shallowSeqDataToAppend) {
            String outputStringForFile =
                    dataToAppend.sampleBarcode() + DELIMITER + dataToAppend.sampleId() + DELIMITER + dataToAppend.purityShallowSeq()
                            + DELIMITER + dataToAppend.hasReliableQuality() + DELIMITER + dataToAppend.hasReliablePurity() + "\n";
            writer.write(outputStringForFile);
            LOGGER.info("Sample barcode: {} is added to shallow seq db!", dataToAppend.sampleBarcode());

        }
        writer.close();
    }

    private static boolean checkInputs(@NotNull CommandLine cmd) {
        String runsDirectory = cmd.getOptionValue(RUNS_DIRECTORY);

        boolean allParamsPresent = !Utils.anyNull(runsDirectory,
                cmd.getOptionValue(SHALLOW_SEQ_TSV),
                cmd.getOptionValue(PURPLE_PURITY_P4_TSV),
                cmd.getOptionValue(PURPLE_PURITY_P5_TSV),
                cmd.getOptionValue(PURPLE_QC_FILE),
                cmd.getOptionValue(PIPELINE_VERSION_FILE));

        boolean validRunDirectories = true;
        if (allParamsPresent) {
            File runDirectoryDb = new File(runsDirectory);

            if (!runDirectoryDb.exists() || !runDirectoryDb.isDirectory()) {
                validRunDirectories = false;
                LOGGER.warn("Shallow seq dir {} does not exist or is not a directory", runDirectoryDb);
            }
        }

        return validRunDirectories && allParamsPresent;
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();

        options.addOption(RUNS_DIRECTORY, true, "Path towards the folder containing all shallow seq runs .");

        options.addOption(SHALLOW_SEQ_TSV, true, "Path towards output file for the shallow seq db TSV.");

        options.addOption(PURPLE_PURITY_P4_TSV, true, "Path towards the purple purity TSV of P4 and lower.");
        options.addOption(PURPLE_PURITY_P5_TSV, true, "Path towards the purple purity TSV of P5.");
        options.addOption(PURPLE_QC_FILE, true, "Path towards the purple qc file.");

        options.addOption(PIPELINE_VERSION_FILE, true, "Path towards the pipeline version");

        return options;
    }
}

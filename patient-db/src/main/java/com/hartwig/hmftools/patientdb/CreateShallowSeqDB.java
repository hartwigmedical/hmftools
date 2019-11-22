package com.hartwig.hmftools.patientdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.ImmutableLimsShallowSeqData;
import com.hartwig.hmftools.common.lims.LimsShallowSeqData;
import com.hartwig.hmftools.common.purple.CheckPurpleQuality;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.qc.PurpleQCFile;
import com.hartwig.hmftools.common.utils.io.reader.LineReader;
import com.hartwig.hmftools.patientdb.context.RunContext;
import com.hartwig.hmftools.patientdb.readers.RunsFolderReader;

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

    private static final String RUNS_DIRECTORY = "runs_dir";

    private static final String SHALLOW_SEQ_CSV = "shallow_seq_csv";
    private static final String PURPLE_PURITY_P4_TSV = "purple_purity_p4_tsv";
    private static final String PURPLE_PURITY_P5_TSV = "purple_purity_p5_tsv";
    private static final String PURPLE_QC_FILE = "purple_qc_file";
    private static final String PIPELINE_VERSION = "pipeline_version_file";

    private static final String PURPLE_DIR = "purple";
    private static final String DELIMITER = ",";

    public static void main(@NotNull final String[] args) throws ParseException, IOException {
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if (!checkInputs(cmd)) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("shallow seq db ", options);
            System.exit(1);
        }

        LOGGER.info("Loading shallow seq runs from {}", cmd.getOptionValue(RUNS_DIRECTORY));
        final List<RunContext> runContexts = loadRunContexts(cmd.getOptionValue(RUNS_DIRECTORY));

        List<LimsShallowSeqData> appendShallowSeqData = extractPurpleFromSample(runContexts,
                cmd.getOptionValue(SHALLOW_SEQ_CSV),
                cmd.getOptionValue(RUNS_DIRECTORY),
                cmd.getOptionValue(PURPLE_QC_FILE),
                cmd.getOptionValue(PURPLE_PURITY_P4_TSV),
                cmd.getOptionValue(PURPLE_PURITY_P5_TSV), cmd.getOptionValue(PIPELINE_VERSION));

        appendToCsv(cmd.getOptionValue(SHALLOW_SEQ_CSV), appendShallowSeqData);

        LOGGER.info("Shallow seq DB is complete!");
    }

    @NotNull
    private static List<LimsShallowSeqData> read(@NotNull String shallowSeqCsv) throws IOException {
        List<String> linesShallowDB = LineReader.build().readLines(new File(shallowSeqCsv).toPath(), line -> line.length() > 0);
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
    private static List<LimsShallowSeqData> extractPurpleFromSample(@NotNull List<RunContext> runContexts,
            @NotNull String shallowSeqOutputCsv, @NotNull String path, @NotNull String purpleQCFile, @NotNull String purplePurityP4,
            @NotNull String purplePurityP5, @NotNull String pipelineVersion) throws IOException {

        List<LimsShallowSeqData> shallowSeqData = read(shallowSeqOutputCsv);

        List<LimsShallowSeqData> appendShallowSeqData = Lists.newArrayList();

        for (RunContext runInfo : runContexts) {
            String tumorSample = runInfo.tumorSample();
            String setPath = path + File.separator + runInfo.setName();
            String sampleBarcode = runInfo.tumorBarcodeSample();

            String purplePurityTsvExt;
            if (new File(setPath + File.separator + pipelineVersion).exists()) {
                purplePurityTsvExt = purplePurityP5;
            } else {
                purplePurityTsvExt = purplePurityP4;
            }
            String purplePurityTsv = setPath + File.separator + PURPLE_DIR + File.separator + tumorSample + purplePurityTsvExt;
            String purpleQcFile = setPath + File.separator + PURPLE_DIR + File.separator + tumorSample + purpleQCFile;

            PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
            PurpleQC purpleQC = PurpleQCFile.read(purpleQcFile);

            boolean hasReliableQuality = CheckPurpleQuality.checkHasReliableQuality(purpleQC);
            boolean hasReliablePurity = CheckPurpleQuality.checkHasReliablePurity(purityContext);
            DecimalFormat decimalFormat = new DecimalFormat("#.##");
            String purity = decimalFormat.format(purityContext.bestFit().purity());

            boolean inFile = false;
            if (shallowSeqData.size() == 0) {
                appendShallowSeqData.add(ImmutableLimsShallowSeqData.builder()
                        .sampleBarcode(sampleBarcode)
                        .sampleId(tumorSample)
                        .purityShallowSeq(purity)
                        .hasReliableQuality(hasReliableQuality)
                        .hasReliablePurity(hasReliablePurity)
                        .build());
                LOGGER.info("Set: {} is added to shallow list!", setPath);
            } else {
                for (LimsShallowSeqData sample : shallowSeqData) {
                    if (sample.sampleBarcode().equals(sampleBarcode)) {
                        LOGGER.warn("Sample barcode is already present in file. Skipping set: {} with sample barcode: {} for"
                                + " writing to shallow seq db!", setPath, sampleBarcode);
                        inFile = true;
                    }
                }
                if (!inFile && !sampleBarcode.equals(Strings.EMPTY)) {
                    appendShallowSeqData.add(ImmutableLimsShallowSeqData.builder()
                            .sampleBarcode(sampleBarcode)
                            .sampleId(tumorSample)
                            .purityShallowSeq(purity)
                            .hasReliableQuality(hasReliableQuality)
                            .hasReliablePurity(hasReliablePurity)
                            .build());
                    LOGGER.info("Set: {} is added to shallow list!", setPath);
                }
            }
        }
        return appendShallowSeqData;
    }

    @NotNull
    private static List<RunContext> loadRunContexts(@NotNull String runsDirectory) throws IOException {
        final List<RunContext> runContexts = RunsFolderReader.extractRunContexts(new File(runsDirectory));
        LOGGER.info(" Loaded run contexts from {} ({} sets)", runsDirectory, runContexts.size());

        return runContexts;
    }

    private static void appendToCsv(@NotNull String shallowSeqCsv, @NotNull List<LimsShallowSeqData> appendingShallowSeqData)
            throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(shallowSeqCsv, true));
        for (LimsShallowSeqData stringToAppend : appendingShallowSeqData) {
            String outputStringForFile =
                    stringToAppend.sampleBarcode() + DELIMITER + stringToAppend.sampleId() + DELIMITER + stringToAppend.purityShallowSeq()
                            + DELIMITER + stringToAppend.hasReliableQuality() + DELIMITER + stringToAppend.hasReliablePurity() + "\n";
            writer.write(outputStringForFile);
            LOGGER.info("Sample barcode: {} is added to shallow seq db!", stringToAppend.sampleBarcode());

        }
        writer.close();
    }

    private static boolean checkInputs(@NotNull CommandLine cmd) {
        final String runsDirectory = cmd.getOptionValue(RUNS_DIRECTORY);

        boolean allParamsPresent = !Utils.anyNull(runsDirectory,
                cmd.getOptionValue(SHALLOW_SEQ_CSV),
                cmd.getOptionValue(PURPLE_PURITY_P4_TSV),
                cmd.getOptionValue(PURPLE_PURITY_P5_TSV),
                cmd.getOptionValue(PURPLE_QC_FILE),
                cmd.getOptionValue(PIPELINE_VERSION));

        boolean validRunDirectories = true;
        if (allParamsPresent) {
            final File runDirectoryDb = new File(runsDirectory);

            if (!runDirectoryDb.exists() || !runDirectoryDb.isDirectory()) {
                validRunDirectories = false;
                LOGGER.warn("Shallow seq dir {} does not exist or is not a directory", runDirectoryDb);
            }
        }

        return validRunDirectories && allParamsPresent;
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

        options.addOption(SHALLOW_SEQ_CSV, true, "Path towards output file for the shallow seq db CSV.");

        options.addOption(PURPLE_PURITY_P4_TSV, true, "Path towards the purple purity TSV of P4 and lower.");
        options.addOption(PURPLE_PURITY_P5_TSV, true, "Path towards the purple purity TSV of P5.");
        options.addOption(PURPLE_QC_FILE, true, "Path towards the purple qc file.");

        options.addOption(PIPELINE_VERSION, true, "Path towards the pipeline version");

        return options;
    }
}

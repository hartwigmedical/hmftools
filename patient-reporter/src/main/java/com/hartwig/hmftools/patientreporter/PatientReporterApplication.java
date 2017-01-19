package com.hartwig.hmftools.patientreporter;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.copynumber.cnv.CNVFileLoader;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.VariantConsequence;
import com.hartwig.hmftools.common.variant.vcf.VCFFileLoader;
import com.hartwig.hmftools.common.variant.vcf.VCFFileWriter;
import com.hartwig.hmftools.common.variant.vcf.VCFSomaticFile;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyser;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberStats;
import com.hartwig.hmftools.patientreporter.slicing.GenomeRegion;
import com.hartwig.hmftools.patientreporter.slicing.Slicer;
import com.hartwig.hmftools.patientreporter.slicing.SlicerFactory;
import com.hartwig.hmftools.patientreporter.util.ConsequenceCount;
import com.hartwig.hmftools.patientreporter.variants.VariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.VariantInterpreter;

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
import org.jetbrains.annotations.Nullable;

public class PatientReporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(PatientReporterApplication.class);
    private static final String SOMATIC_EXTENSION = "_melted.vcf";
    private static final String COPYNUMBER_DIRECTORY = "copyNumber";
    private static final String COPYNUMBER_EXTENSION = ".bam_CNVs";
    private static final String FREEC_DIRECTORY = "freec";

    private static final String RUN_DIRECTORY_ARGS_DESC = "A path towards a single rundir.";
    private static final String RUN_DIRECTORY = "rundir";

    private static final String CPCT_SLICING_BED_ARGS_DESC = "A path towards the CPCT slicing bed.";
    private static final String CPCT_SLICING_BED = "cpct_slicing_bed";

    private static final String HIGH_CONFIDENCE_BED_ARGS_DESC = "A path towards the high confidence bed.";
    private static final String HIGH_CONFIDENCE_BED = "high_confidence_bed";

    private static final String HMF_SLICING_BED_ARGS_DESC = "A path towards the HMF slicing bed.";
    private static final String HMF_SLICING_BED = "hmf_slicing_bed";

    private static final String OUTPUT_DIR_ARGS_DESC = "A path where, if provided, output files will be written to.";
    private static final String OUTPUT_DIR = "output_dir";

    private static final String BATCH_MODE_ARGS_DESC = "If set, runs in batch mode (Caution!!! Korneel Only)";
    private static final String BATCH_MODE = "batch_mode";

    @NotNull
    private final String runDirectory;
    @NotNull
    private final VariantInterpreter variantInterpreter;
    @NotNull
    private final Slicer hmfSlicer;
    @Nullable
    private final String outputDirectory;
    private final boolean batchMode;

    public static void main(final String... args) throws ParseException, IOException, HartwigException {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final String runDir = cmd.getOptionValue(RUN_DIRECTORY);
        final String cpctSlicingBed = cmd.getOptionValue(CPCT_SLICING_BED);
        final String highConfidenceBed = cmd.getOptionValue(HIGH_CONFIDENCE_BED);
        final String hmfSlicingBed = cmd.getOptionValue(HMF_SLICING_BED);
        final String outputDirectory = cmd.getOptionValue(OUTPUT_DIR);
        final boolean batchMode = cmd.hasOption(BATCH_MODE);

        if (runDir == null || cpctSlicingBed == null || highConfidenceBed == null || hmfSlicingBed == null) {
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("Patient-Reporter", options);
            System.exit(1);
        }

        if (outputDirectory != null) {
            final Path outputPath = new File(outputDirectory).toPath();
            if (!Files.exists(outputPath) || !Files.isDirectory(outputPath)) {
                LOGGER.warn(OUTPUT_DIR + " has to be an existing directory!");
                System.exit(1);
            }
        }

        final Slicer hmfSlicingRegion = SlicerFactory.fromBedFile(hmfSlicingBed);
        final VariantInterpreter variantInterpreter = VariantInterpreter.fromSlicingRegions(hmfSlicingRegion,
                SlicerFactory.fromBedFile(highConfidenceBed), SlicerFactory.fromBedFile(cpctSlicingBed));
        new PatientReporterApplication(runDir, variantInterpreter, hmfSlicingRegion, outputDirectory, batchMode).run();
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();

        options.addOption(RUN_DIRECTORY, true, RUN_DIRECTORY_ARGS_DESC);
        options.addOption(CPCT_SLICING_BED, true, CPCT_SLICING_BED_ARGS_DESC);
        options.addOption(HIGH_CONFIDENCE_BED, true, HIGH_CONFIDENCE_BED_ARGS_DESC);
        options.addOption(HMF_SLICING_BED, true, HMF_SLICING_BED_ARGS_DESC);
        options.addOption(OUTPUT_DIR, true, OUTPUT_DIR_ARGS_DESC);
        options.addOption(BATCH_MODE, false, BATCH_MODE_ARGS_DESC);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args)
            throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    PatientReporterApplication(@NotNull final String runDirectory,
            @NotNull final VariantInterpreter variantInterpreter, @NotNull final Slicer hmfSlicer,
            @Nullable final String outputDirectory, final boolean batchMode) {
        this.runDirectory = runDirectory;
        this.variantInterpreter = variantInterpreter;
        this.hmfSlicer = hmfSlicer;
        this.outputDirectory = outputDirectory;
        this.batchMode = batchMode;
    }

    void run() throws IOException, HartwigException {
        if (batchMode) {
            batchRun();
        } else {
            patientRun();
        }
    }

    private void batchRun() throws IOException, HartwigException {
        // KODU: We assume "run directory" is a path with a lot of directories on which we can all run in patient mode.
        final VariantConsequence[] consequences = VariantConsequence.values();

        String header = "SAMPLE,VARIANT_COUNT,PASS_ONLY_COUNT,CONSENSUS_COUNT,MISSENSE_COUNT,CONSEQUENCE_COUNT";
        for (final VariantConsequence consequence : consequences) {
            header += ("," + consequence.name() + "_COUNT");
        }
        System.out.println(header);

        for (final Path run : Files.list(new File(runDirectory).toPath()).collect(Collectors.toList())) {
            final VCFSomaticFile variantFile = loadVariantFile(run.toFile().getPath());
            final VariantAnalysis analysis = variantInterpreter.run(variantFile.variants());

            final Map<VariantConsequence, Integer> counts = ConsequenceCount.count(analysis.consensusPassedVariants());
            String consequenceList = Strings.EMPTY;
            for (final VariantConsequence consequence : consequences) {
                consequenceList += ("," + counts.get(consequence).toString());
            }

            final String out =
                    variantFile.sample() + "," + variantFile.variants().size() + "," + analysis.passedVariants().size()
                            + "," + analysis.consensusPassedVariants().size() + ","
                    + analysis.missenseVariants().size() + "," + analysis.consequencePassedVariants().size()
                    + consequenceList;
            System.out.println(out);
        }
    }

    private void patientRun() throws IOException, HartwigException {
        LOGGER.info("Running patient reporter on " + runDirectory);

        final VCFSomaticFile variantFile = loadVariantFile(runDirectory);
        LOGGER.info("  Total number of variants : " + variantFile.variants().size());

        final VariantAnalysis analysis = variantInterpreter.run(variantFile.variants());
        analyzeCopyNumbers(variantFile.sample());

        LOGGER.info("  Number of variants after applying pass-only filter : " + analysis.passedVariants().size());
        LOGGER.info(
                "  Number of variants after applying consensus rule : " + analysis.consensusPassedVariants().size());
        LOGGER.info("  Number of missense variants in consensus rule (mutational load) : "
                + analysis.missenseVariants().size());
        LOGGER.info("  Number of consequential variants to report : " + analysis.consequencePassedVariants().size());

        if (outputDirectory != null) {
            final String consensusVCF =
                    outputDirectory + File.separator + variantFile.sample() + "_consensus_variants.vcf";
            VCFFileWriter.writeSomaticVCF(consensusVCF, analysis.consensusPassedVariants());
            LOGGER.info("    Written consensus-passed variants to " + consensusVCF);

            final String consequenceVCF =
                    outputDirectory + File.separator + variantFile.sample() + "_consequential_variants.vcf";
            VCFFileWriter.writeSomaticVCF(consequenceVCF, analysis.consequencePassedVariants());
            LOGGER.info("    Written consequential variants to " + consequenceVCF);
        }
    }

    private void analyzeCopyNumbers(@NotNull final String sample) throws IOException, HartwigException {
        final String cnvBasePath = guessCNVBasePath(sample) + File.separator + FREEC_DIRECTORY;
        List<CopyNumber> copyNumbers;

        try {
            copyNumbers = CNVFileLoader.loadCNV(cnvBasePath, sample, COPYNUMBER_EXTENSION);
        } catch (EmptyFileException e) {
            // KODU: It could be that the sample simply does not have any amplifications...
            copyNumbers = Lists.newArrayList();
        }

        final Map<GenomeRegion, CopyNumberStats> stats = CopyNumberAnalyser.run(hmfSlicer.regions(), copyNumbers);
        LOGGER.info("  Determined copy number stats for " + stats.size() + " genomic regions");

        if (outputDirectory != null) {
            final List<String> lines = Lists.newArrayList();
            lines.add("GENE,CNV_MIN,CNV_MEAN,CNV_MAX");
            for (final Map.Entry<GenomeRegion, CopyNumberStats> entry : stats.entrySet()) {
                final CopyNumberStats stat = entry.getValue();
                if (stat.min() != 2 || stat.max() != 2) {
                    lines.add(entry.getKey().annotation() + "," + stat.min() + "," + stat.mean() + "," + stat.max());
                }
            }
            final String filePath = outputDirectory + File.separator + sample + "_CNV.csv";
            Files.write(new File(filePath).toPath(), lines);
            LOGGER.info("    Written all non-default CNV stats to " + filePath);
        }
    }

    @NotNull
    private String guessCNVBasePath(@NotNull final String sample) throws IOException {
        final String basePath = runDirectory + File.separator + COPYNUMBER_DIRECTORY;

        for (Path path : Files.list(new File(basePath).toPath()).collect(Collectors.toList())) {
            if (path.toFile().isDirectory() && path.getFileName().toFile().getName().contains(sample)) {
                return path.toString();
            }
        }

        throw new FileNotFoundException(
                "Could not determine CNV location in " + runDirectory + " using sample " + sample);
    }

    @NotNull
    private static VCFSomaticFile loadVariantFile(@NotNull final String path) throws IOException, HartwigException {
        return VCFFileLoader.loadSomaticVCF(path, SOMATIC_EXTENSION);
    }
}

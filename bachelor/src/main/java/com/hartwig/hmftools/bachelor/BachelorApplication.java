package com.hartwig.hmftools.bachelor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.gene.GeneCopyNumberFile;

import nl.hartwigmedicalfoundation.bachelor.Program;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class BachelorApplication {

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);
    private static final String CONFIG_XML = "configXml";
    private static final String CONFIG_DIRECTORY = "configDirectory";
    private static final String RUN_DIRECTORY = "runDirectory";
    private static final String BATCH_DIRECTORY = "batchDirectory";
    private static final String OUTPUT = "output";
    private static final String VALIDATE = "validate";
    private static final String GERMLINE = "germline";
    private static final String SOMATIC = "somatic";

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(CONFIG_DIRECTORY).required(false).hasArg().desc("folder to find program XMLs").build());
        options.addOption(Option.builder(CONFIG_XML).required(false).hasArg().desc("single config XML to run").build());
        options.addOption(Option.builder(OUTPUT).required().hasArg().desc("output file").build());
        options.addOption(Option.builder(RUN_DIRECTORY).required(false).hasArg().desc("the run directory to look for inputs").build());
        options.addOption(Option.builder(BATCH_DIRECTORY).required(false).hasArg().desc("runs directory to batch process").build());
        options.addOption(Option.builder(VALIDATE).required(false).desc("only validate the configs").build());
        options.addOption(Option.builder(GERMLINE).required(false).desc("process the germline file only").build());
        options.addOption(Option.builder(SOMATIC).required(false).desc("process the somatic file only").build());
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Bachelor", "Determines your eligibility!", options, "", true);
        System.exit(1);
    }

    private static Collection<EligibilityReport> processVCF(final String patient, final boolean isGermline, final File vcf,
            final BachelorEligibility eligibility) {

        final String tag = isGermline ? "germline" : "somatic";
        LOGGER.info("process {} vcf: {}", tag, vcf.getPath());

        try (final VCFFileReader reader = new VCFFileReader(vcf, true)) {
            // TODO: always correct? germline has R,T somatic has just T
            final String sample = reader.getFileHeader().getGenotypeSamples().get(0);
            return eligibility.processVCF(patient, sample, tag, reader);
        } catch (final TribbleException e) {
            LOGGER.error("error with VCF file {}: {}", vcf.getPath(), e.getMessage());
            return Collections.emptyList();
        }
    }

    private static Collection<EligibilityReport> processPurpleCNV(final String patient, final File cnv,
            final BachelorEligibility eligibility) {
        LOGGER.info("process cnv: {}", cnv.getPath());
        try {
            final List<GeneCopyNumber> copyNumbers = GeneCopyNumberFile.read(cnv);
            return eligibility.processCopyNumbers(patient, copyNumbers);
        } catch (final IOException e) {
            LOGGER.error("error with CNV file {}: {}", cnv.getPath(), e.getMessage());
            return Collections.emptyList();
        }
    }

    private static Collection<EligibilityReport> processSV(final File vcf, final BachelorEligibility eligibility) {
        LOGGER.info("process sv: {}", vcf.getPath());
        return Collections.emptyList();
    }

    private static File process(final BachelorEligibility eligibility, final RunDirectory run, final boolean germline,
            final boolean somatic) {

        final String patient = run.getPatientID();
        final boolean doGermline = run.germline != null && germline;
        final boolean doSomatic = run.somatic != null && somatic;

        LOGGER.info("processing run: {}", patient);

        final List<EligibilityReport> result = Lists.newArrayList();
        if (doGermline) {
            result.addAll(processVCF(patient, true, run.germline, eligibility));
        }
        if (doSomatic) {
            result.addAll(processVCF(patient, false, run.somatic, eligibility));
        }
        if (run.copyNumber != null) {
            result.addAll(processPurpleCNV(patient, run.copyNumber, eligibility));
        }
        if (run.structuralVariants != null) {
            result.addAll(processSV(run.structuralVariants, eligibility));
        }

        try {
            final File file = File.createTempFile(patient, ".csv");
            file.deleteOnExit();

            outputToFile(result, file);

            return file;
        } catch (final IOException e) {
            LOGGER.error("error with temporary file for {}", run.prefix);
            return null;
        }
    }

    private static String fileHeader() {
        return String.join(",", Arrays.asList("PATIENT", "SOURCE", "PROGRAM", "ID", "GENES", "CHROM", "POS", "REF", "ALTS", "EFFECTS"));
    }

    private static void outputToFile(final Collection<EligibilityReport> reports, final File outputFile) throws IOException {
        try (final BufferedWriter writer = Files.newBufferedWriter(outputFile.toPath())) {
            for (final EligibilityReport report : reports) {
                for (final VariantModel model : report.variants()) {
                    final VariantContext v = model.Context;

                    final String alts = v.getAlternateAlleles().stream().map(Object::toString).collect(Collectors.joining("|"));
                    final String effects =
                            String.join("|", model.Annotations.stream().flatMap(a -> a.Effects.stream()).collect(Collectors.toSet()));
                    final String genes = String.join("|", model.Annotations.stream().map(a -> a.GeneName).collect(Collectors.toSet()));

                    writer.write(String.format("%s,%s,%s,%s,%s,%s,%d,%s,%s,%s", report.patient(), report.tag(), report.program(), v.getID(),
                            genes, v.getContig(), v.getStart(), v.getReference(), alts, effects));
                    writer.newLine();
                }
            }
            writer.close();
        }
    }

    public static void main(final String... args) {
        final Options options = createOptions();
        try {
            final CommandLine cmd = createCommandLine(options, args);

            // load configs
            final Map<String, Program> map;
            if (cmd.hasOption(CONFIG_DIRECTORY)) {
                map = BachelorHelper.loadXML(Paths.get(cmd.getOptionValue(CONFIG_DIRECTORY)));
            } else if (cmd.hasOption(CONFIG_XML)) {
                map = BachelorHelper.loadXML(Paths.get(cmd.getOptionValue(CONFIG_XML)));
            } else {
                LOGGER.error("config directory or xml required!");
                System.exit(1);
                return;
            }

            if (cmd.hasOption(VALIDATE)) {
                System.exit(0);
                return;
            }
            if (map.isEmpty()) {
                LOGGER.error("no programs loaded!");
                System.exit(1);
                return;
            }

            final BachelorEligibility eligibility = BachelorEligibility.fromMap(map);

            LOGGER.info("beginning processing...");

            final boolean germline = !cmd.hasOption(SOMATIC);
            final boolean somatic = !cmd.hasOption(GERMLINE);
            if (!(germline || somatic)) {
                LOGGER.error("can't have -germline and -somatic set at the same time");
            }

            final List<File> filesToMerge;
            if (cmd.hasOption(BATCH_DIRECTORY)) {
                final Path root = Paths.get(cmd.getOptionValue(BATCH_DIRECTORY));
                try (final Stream<Path> stream = Files.walk(root, 1, FileVisitOption.FOLLOW_LINKS).parallel()) {
                    filesToMerge = stream.filter(p -> p.toFile().isDirectory())
                            .filter(p -> !p.equals(root))
                            .map(RunDirectory::new)
                            .map(run -> process(eligibility, run, germline, somatic))
                            .collect(Collectors.toList());
                }
            } else if (cmd.hasOption(RUN_DIRECTORY)) {
                filesToMerge = Collections.singletonList(
                        process(eligibility, new RunDirectory(Paths.get(cmd.getOptionValue(RUN_DIRECTORY))), germline, somatic));
            } else {
                LOGGER.error("requires either a batch or single run directory");
                System.exit(1);
                return;
            }

            LOGGER.info("... processing complete!");
            LOGGER.info("merging to CSV {} ...", cmd.getOptionValue(OUTPUT));

            // TODO: better way to join files? using FileChannels?

            try (final BufferedWriter writer = Files.newBufferedWriter(Paths.get(cmd.getOptionValue(OUTPUT)))) {

                // header
                writer.write(fileHeader());
                writer.newLine();

                for (final File file : filesToMerge) {
                    final List<String> lines = Files.readAllLines(file.toPath());
                    for (final String line : lines) {
                        writer.write(line);
                        writer.newLine();
                    }
                }

            }

            LOGGER.info("... output complete!");

            LOGGER.info("bachelor done!");

        } catch (final ParseException e) {
            printHelpAndExit(options);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

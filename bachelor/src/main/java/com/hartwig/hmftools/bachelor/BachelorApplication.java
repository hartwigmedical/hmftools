package com.hartwig.hmftools.bachelor;

import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;

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

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class BachelorApplication {

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);
    private static final String CONFIG_DIRECTORY = "configDirectory";
    private static final String RUN_DIRECTORY = "runDirectory";
    private static final String BATCH_DIRECTORY = "batchDirectory";
    private static final String OUTPUT = "output";

    // DEBUG
    private static final String GERMLINE_VCF = "germlineVcf";
    private static final String SOMATIC_VCF = "somaticVcf";
    private static final String PURPLE_FILE = "purple";
    private static final String BPI_VCF = "bpiVcf";

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(CONFIG_DIRECTORY).required().hasArg().desc("folder to find program XMLs").build());
        options.addOption(Option.builder(RUN_DIRECTORY).required(false).hasArg().desc("the run directory to look for inputs").build());
        options.addOption(Option.builder(OUTPUT).hasArg().desc("output file").build());
        options.addOption(Option.builder(GERMLINE_VCF).required(false).hasArg().desc("germline vcf to process").build());
        options.addOption(Option.builder(SOMATIC_VCF).required(false).hasArg().desc("somatic vcf to process").build());
        options.addOption(Option.builder(BPI_VCF).required(false).hasArg().desc("BPI structural variant vcf to process").build());
        options.addOption(Option.builder(PURPLE_FILE).required(false).hasArg().desc("purple cnv file to process").build());
        options.addOption(Option.builder(BATCH_DIRECTORY).required(false).hasArg().desc("runs directory to batch process").build());
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

    private static void processVCF(final boolean isGermline, final File vcf, final BachelorEligibility eligibility,
            final List<EligibilityReport> reports) throws Exception {
        LOGGER.info("process vcf: {}", vcf.getPath());

        final VCFFileReader reader = new VCFFileReader(vcf, false);

        // TODO: always correct? germline has R,T somatic has just T
        final String sample = reader.getFileHeader().getGenotypeSamples().get(0);
        final Collection<EligibilityReport> result = eligibility.processVCF(sample, isGermline ? "germline" : "somatic", reader);
        reports.addAll(result);

        reader.close();
    }

    private static void processPurpleCNV(final File cnv, final BachelorEligibility eligibility, final List<EligibilityReport> reports) {
        LOGGER.info("process cnv: {}", cnv.getPath());
    }

    private static void processSV(final File vcf, final BachelorEligibility eligibility, final List<EligibilityReport> reports) {
        LOGGER.info("process sv: {}", vcf.getPath());
    }

    private static void process(final BachelorEligibility eligibility, final List<EligibilityReport> reports, final File germline,
            final File somatic, final File copyNumber, final File structuralVariants) throws Exception {
        if (germline != null) {
            processVCF(true, germline, eligibility, reports);
        }
        if (somatic != null) {
            processVCF(false, somatic, eligibility, reports);
        }
        if (copyNumber != null) {
            processPurpleCNV(copyNumber, eligibility, reports);
        }
        if (structuralVariants != null) {
            processSV(structuralVariants, eligibility, reports);
        }
    }

    public static void main(final String... args) {
        final Options options = createOptions();
        try {
            final CommandLine cmd = createCommandLine(options, args);

            final Path configPath = Paths.get(cmd.getOptionValue(CONFIG_DIRECTORY));
            final Map<String, Program> map = BachelorHelper.loadXML(configPath);
            final BachelorEligibility eligibility = BachelorEligibility.fromMap(map);
            final List<EligibilityReport> reports = Lists.newArrayList();

            final List<RunDirectory> runDirectories = Lists.newArrayList();
            if (cmd.hasOption(BATCH_DIRECTORY)) {
                final File[] runs = Paths.get(cmd.getOptionValue(BATCH_DIRECTORY)).toFile().listFiles(File::isDirectory);
                if (runs != null) {
                    Stream.of(runs).forEach(r -> runDirectories.add(new RunDirectory(r.toPath())));
                }
            } else if (cmd.hasOption(RUN_DIRECTORY)) {
                runDirectories.add(new RunDirectory(Paths.get(cmd.getOptionValue(RUN_DIRECTORY))));
            } else {
                process(eligibility, reports, cmd.hasOption(GERMLINE_VCF) ? Paths.get(cmd.getOptionValue(GERMLINE_VCF)).toFile() : null,
                        cmd.hasOption(SOMATIC_VCF) ? Paths.get(cmd.getOptionValue(SOMATIC_VCF)).toFile() : null,
                        cmd.hasOption(PURPLE_FILE) ? Paths.get(cmd.getOptionValue(PURPLE_FILE)).toFile() : null,
                        cmd.hasOption(BPI_VCF) ? Paths.get(cmd.getOptionValue(BPI_VCF)).toFile() : null);
            }

            for (final RunDirectory run : runDirectories) {
                LOGGER.info("processing run: {}", run.prefix);

                final File germline = run.findGermline();
                final File somatic = run.findSomatic();
                final File copyNumber = run.findCopyNumber();
                final File structuralVariants = run.findStructuralVariants();

                process(eligibility, reports, germline, somatic, copyNumber, structuralVariants);
            }

            // output results

            if (cmd.hasOption(OUTPUT)) {
                final BufferedWriter writer = Files.newBufferedWriter(Paths.get(cmd.getOptionValue(OUTPUT)));

                // header
                writer.write(String.join(",", Arrays.asList("SAMPLE", "SOURCE", "PROGRAM", "ID", "CHROM", "POS", "REF", "ALTS")));
                writer.newLine();

                // data
                for (final EligibilityReport report : reports) {
                    for (final VariantContext v : report.variants()) {
                        final String alts =
                                String.join("|", v.getAlternateAlleles().stream().map(Object::toString).collect(Collectors.toList()));
                        writer.write(String.format("%s,%s,%s,%s,%s,%d,%s,%s", report.sample(), report.tag(), report.program(), v.getID(),
                                v.getContig(), v.getStart(), v.getReference(), alts));
                        writer.newLine();
                    }
                }

                writer.close();
            }

        } catch (final ParseException e) {
            printHelpAndExit(options);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

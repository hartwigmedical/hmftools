package com.hartwig.hmftools.bachelor;

import java.io.BufferedWriter;
import java.io.File;
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
    private static final String CONFIG_DIRECTORY = "configDirectory";
    private static final String RUN_DIRECTORY = "runDirectory";
    private static final String BATCH_DIRECTORY = "batchDirectory";
    private static final String OUTPUT = "output";
    private static final String VALIDATE = "validate";

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(CONFIG_DIRECTORY).required().hasArg().desc("folder to find program XMLs").build());
        options.addOption(Option.builder(OUTPUT).required().hasArg().desc("output file").build());
        options.addOption(Option.builder(RUN_DIRECTORY).required(false).hasArg().desc("the run directory to look for inputs").build());
        options.addOption(Option.builder(BATCH_DIRECTORY).required(false).hasArg().desc("runs directory to batch process").build());
        options.addOption(Option.builder(VALIDATE).required(false).desc("only validate the configs").build());
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

    private static Collection<EligibilityReport> processVCF(final boolean isGermline, final File vcf,
            final BachelorEligibility eligibility) {
        LOGGER.info("process {} vcf: {}", isGermline ? "gemline" : "somatic", vcf.getPath());
        try (final VCFFileReader reader = new VCFFileReader(vcf, false)) {
            // TODO: always correct? germline has R,T somatic has just T
            final String sample = reader.getFileHeader().getGenotypeSamples().get(0);
            return eligibility.processVCF(sample, isGermline ? "germline" : "somatic", reader);
        } catch (final TribbleException e) {
            LOGGER.error("error with VCF file {}: {}", vcf.getPath(), e.getMessage());
            return Collections.emptyList();
        }
    }

    private static Collection<EligibilityReport> processPurpleCNV(final File cnv, final BachelorEligibility eligibility) {
        LOGGER.info("process cnv: {}", cnv.getPath());
        return Collections.emptyList();
    }

    private static Collection<EligibilityReport> processSV(final File vcf, final BachelorEligibility eligibility) {
        LOGGER.info("process sv: {}", vcf.getPath());
        return Collections.emptyList();
    }

    private static Collection<EligibilityReport> process(final BachelorEligibility eligibility, final RunDirectory run) {
        LOGGER.info("processing run: {}", run.prefix);
        final List<EligibilityReport> result = Lists.newArrayList();
        if (run.germline != null) {
            result.addAll(processVCF(true, run.germline, eligibility));
        }
        if (run.somatic != null) {
            result.addAll(processVCF(false, run.somatic, eligibility));
        }
        if (run.copyNumber != null) {
            result.addAll(processPurpleCNV(run.copyNumber, eligibility));
        }
        if (run.structuralVariants != null) {
            result.addAll(processSV(run.structuralVariants, eligibility));
        }
        return result;
    }

    public static void main(final String... args) {
        final Options options = createOptions();
        try {
            final CommandLine cmd = createCommandLine(options, args);

            // load configs

            final Path configPath = Paths.get(cmd.getOptionValue(CONFIG_DIRECTORY));
            final Map<String, Program> map = BachelorHelper.loadXML(configPath);
            if (cmd.hasOption(VALIDATE)) {
                System.exit(0);
                return;
            }

            final BachelorEligibility eligibility = BachelorEligibility.fromMap(map);

            final List<RunDirectory> runDirectories;
            if (cmd.hasOption(BATCH_DIRECTORY)) {
                LOGGER.info("collating files for processing...");
                final Path root = Paths.get(cmd.getOptionValue(BATCH_DIRECTORY));
                try (final Stream<Path> stream = Files.walk(root, 1, FileVisitOption.FOLLOW_LINKS).parallel()) {
                    runDirectories = stream.filter(p -> p.toFile().isDirectory())
                            .filter(p -> !p.equals(root))
                            .map(RunDirectory::new)
                            .collect(Collectors.toList());
                }
            } else if (cmd.hasOption(RUN_DIRECTORY)) {
                runDirectories = Collections.singletonList(new RunDirectory(Paths.get(cmd.getOptionValue(RUN_DIRECTORY))));
            } else {
                LOGGER.error("requires either a batch or single run directory");
                System.exit(1);
                return;
            }

            LOGGER.info("beginning processing...");
            final List<EligibilityReport> reports =
                    runDirectories.parallelStream().flatMap(run -> process(eligibility, run).stream()).collect(Collectors.toList());
            LOGGER.info("... processing complete!");

            // output results

            LOGGER.info("outputting to CSV {} ...", cmd.getOptionValue(OUTPUT));

            try (final BufferedWriter writer = Files.newBufferedWriter(Paths.get(cmd.getOptionValue(OUTPUT)))) {

                // header
                writer.write(
                        String.join(",", Arrays.asList("SAMPLE", "SOURCE", "PROGRAM", "ID", "CHROM", "POS", "REF", "ALTS", "EFFECTS")));
                writer.newLine();

                // data
                for (final EligibilityReport report : reports) {
                    for (final VariantModel model : report.variants()) {

                        final VariantContext v = model.Context;
                        final String alts = v.getAlternateAlleles().stream().map(Object::toString).collect(Collectors.joining("|"));
                        final String effects = model.Annotations.stream().flatMap(a -> a.Effects.stream()).collect(Collectors.joining("|"));

                        writer.write(String.format("%s,%s,%s,%s,%s,%d,%s,%s,%s", report.sample(), report.tag(), report.program(), v.getID(),
                                v.getContig(), v.getStart(), v.getReference(), alts, effects));
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

package com.hartwig.hmftools.bachelor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

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
import org.xml.sax.SAXException;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class BachelorApplication {

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);
    private static final String CONFIG_DIRECTORY = "configDirectory";
    private static final String RUN_DIRECTORY = "runDirectory";
    private static final String OUTPUT = "output";

    // DEBUG
    private static final String GERMLINE_VCF = "germlineVcf";
    private static final String SOMATIC_VCF = "somaticVcf";

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(CONFIG_DIRECTORY).required().hasArg().desc("folder to find program XMLs").build());
        options.addOption(Option.builder(RUN_DIRECTORY).required(false).hasArg().desc("the run directory to look for inputs").build());
        options.addOption(Option.builder(OUTPUT).hasArg().desc("output file").build());
        options.addOption(Option.builder(GERMLINE_VCF).required(false).hasArg().desc("germline vcf to process").build());
        options.addOption(Option.builder(SOMATIC_VCF).required(false).hasArg().desc("somatic vcf to process").build());
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

    private static void processVCF(final String tag, final File vcf, final BachelorEligibility eligibility,
            final List<EligibilityReport> reports, final Map<String, Integer> merged) {
        LOGGER.info("process vcf: {}", vcf.getPath());
        final VCFFileReader reader = new VCFFileReader(vcf, false);
        final Collection<EligibilityReport> result = eligibility.processVCF(tag, reader);
        result.forEach(report -> merged.merge(report.program(), report.variants().size(), (a, b) -> a + b));
        reports.addAll(result);
        reader.close();
    }

    public static void main(final String... args) {
        final Options options = createOptions();
        try {
            final CommandLine cmd = createCommandLine(options, args);

            final Path configPath = Paths.get(cmd.getOptionValue(CONFIG_DIRECTORY));
            final Map<String, Program> map = BachelorHelper.loadXML(configPath);
            final BachelorEligibility eligibility = BachelorEligibility.fromMap(map);

            final File germline;
            final File somatic;

            if (cmd.hasOption(RUN_DIRECTORY)) {
                final RunDirectory run = new RunDirectory(Paths.get(cmd.getOptionValue(RUN_DIRECTORY)));
                germline = run.findGermline();
                somatic = run.findSomatic();
            } else {
                germline = cmd.hasOption(GERMLINE_VCF) ? Paths.get(cmd.getOptionValue(GERMLINE_VCF)).toFile() : null;
                somatic = cmd.hasOption(SOMATIC_VCF) ? Paths.get(cmd.getOptionValue(SOMATIC_VCF)).toFile() : null;
            }

            final List<EligibilityReport> reports = Lists.newArrayList();
            final Map<String, Integer> merged = Maps.newHashMap();

            if (germline != null) {
                processVCF("germline", germline, eligibility, reports, merged);
            }
            if (somatic != null) {
                processVCF("somatic", somatic, eligibility, reports, merged);
            }

            // output results

            merged.entrySet()
                    .stream()
                    .sorted(Comparator.comparingInt(Map.Entry::getValue))
                    .forEach(e -> LOGGER.info("{} = {} variants", e.getKey(), e.getValue()));

            if (cmd.hasOption(OUTPUT)) {
                final BufferedWriter writer = Files.newBufferedWriter(Paths.get(cmd.getOptionValue(OUTPUT)));
                for (final EligibilityReport report : reports) {
                    for (final VariantContext v : report.variants()) {
                        writer.write(
                                String.format("%s,%s,%s,%d,%s,%s" + System.lineSeparator(), report.tag(), report.program(), v.getContig(),
                                        v.getStart(), v.getReference(), String.join("|",
                                                v.getAlternateAlleles().stream().map(Object::toString).collect(Collectors.toList()))));
                    }
                }
                writer.close();
            }

        } catch (final ParseException e) {
            printHelpAndExit(options);
        } catch (SAXException | IOException e) {
            e.printStackTrace();
        }
    }
}

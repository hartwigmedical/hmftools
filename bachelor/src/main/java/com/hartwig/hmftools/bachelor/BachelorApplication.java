package com.hartwig.hmftools.bachelor;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

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

import htsjdk.variant.vcf.VCFFileReader;

public class BachelorApplication {

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);
    private static final String CONFIG_DIRECTORY = "configDirectory";
    private static final String VCF = "vcf";

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(CONFIG_DIRECTORY).required().hasArg().desc("folder to find program XMLs").build());
        options.addOption(Option.builder(VCF).required().numberOfArgs(Option.UNLIMITED_VALUES).desc("vcf to process").build());
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

    private static Map<String, Program> loadXML(final Path path) throws IOException, SAXException {
        final BachelorSchema schema = BachelorSchema.make();

        final List<Program> programs = Files.walk(path)
                .filter(p -> p.toString().endsWith(".xml"))
                .map(schema::processXML)
                .filter(Objects::nonNull)
                .collect(Collectors.toList());

        final Map<String, Program> result = Maps.newHashMap();
        for (final Program p : programs) {
            if (result.containsKey(p.getName())) {
                LOGGER.error("duplicate programs detected: {}", p.getName());
                System.exit(1);
            } else {
                result.put(p.getName(), p);
            }
        }

        return result;
    }

    public static void main(final String... args) {
        final Options options = createOptions();
        try {
            final CommandLine cmd = createCommandLine(options, args);

            final Path configPath = Paths.get(cmd.getOptionValue(CONFIG_DIRECTORY));
            final Map<String, Program> map = loadXML(configPath);
            final BachelorEligibility eligibility = BachelorEligibility.fromMap(map);

            final List<File> vcfFiles =
                    Arrays.stream(cmd.getOptionValues(VCF)).map(s -> Paths.get(s).toFile()).collect(Collectors.toList());

            final Map<String, Integer> merged = Maps.newHashMap();
            for (final File vcf : vcfFiles) {
                LOGGER.info("process vcf: {}", vcf.getPath());
                final VCFFileReader reader = new VCFFileReader(vcf, false);
                final Map<String, Integer> result = eligibility.processVCF(reader);
                result.forEach((k, v) -> merged.merge(k, v, (a, b) -> a + b));
                reader.close();
            }
            LOGGER.info(merged.entrySet().stream().sorted(Comparator.comparingInt(Map.Entry::getValue)).collect(Collectors.toList()));
        } catch (final ParseException e) {
            printHelpAndExit(options);
        } catch (SAXException | IOException e) {
            e.printStackTrace();
        }
    }
}

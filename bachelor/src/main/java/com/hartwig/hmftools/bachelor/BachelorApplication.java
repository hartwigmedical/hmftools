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
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;

import nl.hartwigmedicalfoundation.bachelor.Program;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;

public class BachelorApplication {

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);
    private static final String CONFIG_XML = "configXml";
    private static final String CONFIG_DIRECTORY = "configDirectory";
    private static final String RUN_DIRECTORY = "runDirectory";
    private static final String BATCH_DIRECTORY = "batchDirectory";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String VALIDATE = "validate";
    private static final String GERMLINE = "germline";
    private static final String SOMATIC = "somatic";
    private static final String COPYNUMBER = "copyNumber";
    private static final String SV = "structuralVariants";
    private static final String SAMPLE = "sample";
    private static final String LOG_DEBUG = "log_debug";

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(CONFIG_DIRECTORY).required(false).hasArg().desc("folder to find program XMLs").build());
        options.addOption(Option.builder(CONFIG_XML).required(false).hasArg().desc("single config XML to run").build());
        options.addOption(Option.builder(OUTPUT_DIR).required().hasArg().desc("output file").build());
        options.addOption(Option.builder(RUN_DIRECTORY).required(false).hasArg().desc("the run directory to look for inputs").build());
        options.addOption(Option.builder(BATCH_DIRECTORY).required(false).hasArg().desc("runs directory to batch process").build());
        options.addOption(Option.builder(VALIDATE).required(false).desc("only validate the configs").build());
        options.addOption(Option.builder(GERMLINE).required(false).desc("process the germline file").build());
        options.addOption(Option.builder(SOMATIC).required(false).desc("process the somatic file").build());
        options.addOption(Option.builder(COPYNUMBER).required(false).desc("process the copy number file").build());
        options.addOption(Option.builder(SV).required(false).desc("process the sv file").build());
        options.addOption(Option.builder(SAMPLE).required(false).hasArg().desc("sample id").build());
        options.addOption(Option.builder(LOG_DEBUG).required(false).desc("Sets log level to Debug, off by default").build());
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Bachelor", "Determines your eligibility", options, "", true);
        System.exit(1);
    }

    private static Collection<EligibilityReport> processVCF(final String patient, final boolean isGermline, final File vcf,
            final BachelorEligibility eligibility) {

        final EligibilityReport.ReportType type =
                isGermline ? EligibilityReport.ReportType.GERMLINE_MUTATION : EligibilityReport.ReportType.SOMATIC_MUTATION;

        LOGGER.info("processing {} vcf: {}", type, vcf.getPath());

        try (final VCFFileReader reader = new VCFFileReader(vcf, true)) {
            // assume that the first sample is the germline
            final String sample = reader.getFileHeader().getGenotypeSamples().get(0);
            return eligibility.processVCF(patient, sample, type, reader);
        } catch (final TribbleException e) {
            LOGGER.error("error with VCF file {}: {}", vcf.getPath(), e.getMessage());
            return Collections.emptyList();
        }
    }

    private static Collection<EligibilityReport> processPurpleCNV(final String patient, final File cnv,
            final BachelorEligibility eligibility) {
        LOGGER.info("processing cnv: {}", cnv.getPath());
        try {
            final List<GeneCopyNumber> copyNumbers = GeneCopyNumberFile.read(cnv);
            return eligibility.processCopyNumbers(patient, copyNumbers);
        } catch (final IOException e) {
            LOGGER.error("error with CNV file {}: {}", cnv.getPath(), e.getMessage());
            return Collections.emptyList();
        }
    }

    private static Collection<EligibilityReport> processSV(final String patient, final File vcf, final BachelorEligibility eligibility) {
        LOGGER.info("processing sv: {}", vcf.getPath());

        final StructuralVariantFactory factory = new StructuralVariantFactory(true);
        try {
            try (final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(vcf.getPath(),
                    new VCFCodec(),
                    false)) {
                reader.iterator().forEach(factory::addVariantContext);
            }
        } catch (IOException e) {
            LOGGER.error("error with SV file {}: {}", vcf.getPath(), e.getMessage());
            return Collections.emptyList();
        }

        return eligibility.processStructuralVariants(patient, factory.results());
    }

    private static void process(
            final BachelorEligibility eligibility, final RunDirectory run, final String sampleId,
            final boolean germline, final boolean somatic, final boolean copyNumber, final boolean structuralVariants,
            final BufferedWriter allDataWriter, final BufferedWriter bedFileWriter)
    {
        final String patient = run.getPatientID();
        final boolean doGermline = run.germline() != null && germline;
        final boolean doSomatic = run.somatic() != null && somatic;
        final boolean doCopyNumber = run.copyNumber() != null && copyNumber;
        final boolean doStructuralVariants = run.structuralVariants() != null && structuralVariants;

        LOGGER.info("processing run for patient({}) from file({})", patient, run.prefix());

        final List<EligibilityReport> result = Lists.newArrayList();
        if (doGermline) {
            result.addAll(processVCF(patient, true, run.germline(), eligibility));
        }
        if (doSomatic) {
            result.addAll(processVCF(patient, false, run.somatic(), eligibility));
        }
        if (doCopyNumber) {
            result.addAll(processPurpleCNV(patient, run.copyNumber(), eligibility));
        }
        if (doStructuralVariants) {
            result.addAll(processSV(patient, run.structuralVariants(), eligibility));
        }

        try {
            for (final EligibilityReport r : result) {

                allDataWriter.write(String.format("%s,%s,%s,%s,%s,%s,%s,%d,%s,%s,%s,%s",
                        sampleId, r.source().toString(), r.program(),  r.id(),
                        r.genes(), r.transcriptId(), r.chrom(), r.pos(), r.ref(), r.alts(), r.effects(), r.hgvsProtein()));
                allDataWriter.newLine();

                bedFileWriter.write(String.format("%s\t%s\t%d\t%d",
                        sampleId, r.chrom(), r.pos() - 1, r.pos()));
                bedFileWriter.newLine();
            }
        }
        catch (final IOException e) {
            LOGGER.error("error writing output");
            return;
        }
    }

    private static String fileHeader() {
        return String.join(",",
                Arrays.asList("SAMPLEID", "SOURCE", "PROGRAM", "ID", "GENE", "TRANSCRIPT_ID", "CHROM", "POS", "REF", "ALTS", "EFFECTS", "HGVS_PROTEIN"));
    }

    public static void main(final String... args) {
        final Options options = createOptions();

        try
        {
            final CommandLine cmd = createCommandLine(options, args);

            if (cmd.hasOption(LOG_DEBUG)) {
                Configurator.setRootLevel(Level.DEBUG);
            }

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
                LOGGER.error("no programs loaded, exiting");
                System.exit(1);
                return;
            }

            boolean isBatchRun = cmd.hasOption(BATCH_DIRECTORY);
            boolean isSingleRun = cmd.hasOption(RUN_DIRECTORY);
            final String sampleId = cmd.getOptionValue(SAMPLE);

            if (!isBatchRun && !isSingleRun)
            {
                LOGGER.error("requires either a batch or single run directory");
                System.exit(1);
                return;
            }
            else if(isSingleRun && (sampleId == null || sampleId.isEmpty()))
            {
                LOGGER.error("single run requires sample to be specified");
                System.exit(1);
                return;
            }

            final BachelorEligibility eligibility = BachelorEligibility.fromMap(map);

            if(isBatchRun)
                LOGGER.info("beginning batch run");
            else
                LOGGER.info("beginning single sample run: {}", sampleId);

            final boolean germline = cmd.hasOption(GERMLINE);
            final boolean somatic = cmd.hasOption(SOMATIC);
            final boolean copyNumber = cmd.hasOption(COPYNUMBER);
            final boolean structuralVariants = cmd.hasOption(SV);
            final boolean doAll = !(germline || somatic || copyNumber || structuralVariants);

            try {

                String outputDir = cmd.getOptionValue(OUTPUT_DIR);

                if(!outputDir.endsWith("/"))
                    outputDir += "/";

                String mainFileName = outputDir + "bachelor_output.csv";
                String bedFileName = outputDir + "bachelor_bed.csv";

                final BufferedWriter mainDataWriter = Files.newBufferedWriter(Paths.get(mainFileName));
                mainDataWriter.write(fileHeader());
                mainDataWriter.newLine();

                final BufferedWriter bedFileWriter = Files.newBufferedWriter(Paths.get(bedFileName));

                if (cmd.hasOption(BATCH_DIRECTORY))
                {
                    final Path root = Paths.get(cmd.getOptionValue(BATCH_DIRECTORY));

                    try (final Stream<Path> stream = Files.walk(root, 1, FileVisitOption.FOLLOW_LINKS).parallel())
                    {
                        final List<RunDirectory> runDirectories = stream
                                .filter(p -> p.toFile().isDirectory())
                                .filter(p -> !p.equals(root))
                                .map(RunDirectory::new)
                                .collect(Collectors.toList());

                        LOGGER.info("found {} batch directories", runDirectories.size());

                            // add the filtered and passed SV entries for each file
                            for (final RunDirectory runDir: runDirectories)
                            {
                                process(eligibility, runDir, runDir.getPatientID(),
                                        germline || doAll,
                                        somatic || doAll,
                                        copyNumber || doAll,
                                        structuralVariants || doAll,
                                        mainDataWriter, bedFileWriter);
                            }
                    }
                    catch (Exception e)
                    {
                        LOGGER.error("failed walking batch directories");
                    }
                }
                else if (cmd.hasOption(RUN_DIRECTORY))
                {
                    final Path path = Paths.get(cmd.getOptionValue(RUN_DIRECTORY));
                    if (!Files.exists(path)) {
                        LOGGER.error("-runDirectory path does not exist");
                        System.exit(1);
                        return;
                    }
                    process(eligibility, new RunDirectory(path), sampleId,
                            germline || doAll,
                            somatic || doAll,
                            copyNumber || doAll,
                            structuralVariants || doAll,
                            mainDataWriter, bedFileWriter);
                }

                mainDataWriter.close();
                bedFileWriter.close();
            }
            catch(IOException e)
            {
                LOGGER.error("failed writing output");
            }

            LOGGER.info("processing complete");

            LOGGER.info("run complete");

        } catch (final ParseException e) {
            printHelpAndExit(options);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}

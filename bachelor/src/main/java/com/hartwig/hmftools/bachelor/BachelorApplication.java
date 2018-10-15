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
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
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

    private Map<String, Program> mProgramMap;
    private BufferedWriter mMainDataWriter;
    private BufferedWriter mBedFileWriter;

    // config
    private String mOutputDir;
    private String mBatchDirectoy;
    private String mRunDirectoy;
    private boolean mIsBatchRun;
    private boolean mIsSingleRun;
    private String mSampleId;
    private boolean mRunGermline;
    private boolean mRunSomatic;
    private boolean mRunCopyNumber;
    private boolean mRunStructuralVariants;

    public BachelorApplication()
    {
        mProgramMap = null;
        mMainDataWriter = null;
        mBedFileWriter = null;

        mIsSingleRun = false;
        mIsBatchRun = false;
    }

    @NotNull
    private static Options createOptions()
    {
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
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        try
        {
            // load configs
            if (cmd.hasOption(CONFIG_DIRECTORY))
            {
                mProgramMap = BachelorHelper.loadXML(Paths.get(cmd.getOptionValue(CONFIG_DIRECTORY)));
            }
            else if (cmd.hasOption(CONFIG_XML))
            {
                mProgramMap = BachelorHelper.loadXML(Paths.get(cmd.getOptionValue(CONFIG_XML)));
            }
            else
            {
                LOGGER.error("config directory or xml required!");
                return false;
            }
        }
        catch(Exception e)
        {
            LOGGER.error("error loading XML: {}", e.toString());
            return false;
        }

        if (cmd.hasOption(VALIDATE))
        {
            return true;
        }

        if (mProgramMap.isEmpty())
        {
            LOGGER.error("no programs loaded, exiting");
            return false;
        }

        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);
        mSampleId = cmd.getOptionValue(SAMPLE);

        mIsBatchRun = cmd.hasOption(BATCH_DIRECTORY);
        mIsSingleRun = cmd.hasOption(RUN_DIRECTORY);

        if(cmd.hasOption(BATCH_DIRECTORY))
        {
            mBatchDirectoy = cmd.getOptionValue(BATCH_DIRECTORY);
            mIsBatchRun = true;
        }
        else if(cmd.hasOption(RUN_DIRECTORY))
        {
            mIsSingleRun = true;
            mRunDirectoy = cmd.getOptionValue(RUN_DIRECTORY);
        }

        if (!mIsBatchRun && !mIsSingleRun)
        {
            LOGGER.error("requires either a batch or single run directory");
            return false;
        }
        else if (mIsSingleRun && (mSampleId == null || mSampleId.isEmpty()))
        {
            LOGGER.error("single run requires sample to be specified");
            return false;
        }

        mRunGermline = cmd.hasOption(GERMLINE);
        mRunSomatic = cmd.hasOption(SOMATIC);
        mRunCopyNumber = cmd.hasOption(COPYNUMBER);
        mRunStructuralVariants = cmd.hasOption(SV);

        return true;
    }

    public boolean run()
    {
        final BachelorEligibility eligibility = BachelorEligibility.fromMap(mProgramMap);

        if (mIsBatchRun)
        {
            LOGGER.info("beginning batch run");
        }
        else
        {
            LOGGER.info("beginning single sample run: {}", mSampleId);
        }

        try
        {
            if (mIsBatchRun)
            {
                final Path root = Paths.get(mBatchDirectoy);

                try (final Stream<Path> stream = Files.walk(root, 1, FileVisitOption.FOLLOW_LINKS).parallel())
                {
                    final List<RunDirectory> runDirectories = stream.filter(p -> p.toFile().isDirectory())
                            .filter(p -> !p.equals(root))
                            .map(RunDirectory::new)
                            .collect(Collectors.toList());

                    LOGGER.info("found {} batch directories", runDirectories.size());

                    // add the filtered and passed SV entries for each file
                    for (final RunDirectory runDir : runDirectories)
                    {
                        process(eligibility, runDir, runDir.getPatientID());
                    }
                }
                catch (Exception e)
                {
                    LOGGER.error("failed walking batch directories");
                }
            }
            else if (mIsSingleRun)
            {
                final Path path = Paths.get(mRunDirectoy);

                if (!Files.exists(path))
                {
                    LOGGER.error("-runDirectory path does not exist");
                    return false;
                }

                process(eligibility, new RunDirectory(path), mSampleId);
            }

            LOGGER.info("run complete");

            if(mMainDataWriter != null)
                mMainDataWriter.close();

            if(mBedFileWriter != null)
                mBedFileWriter.close();
        }
        catch (IOException e)
        {
            LOGGER.error("failed writing output: {}", e.toString());
        }

        return true;
    }

    private static Collection<EligibilityReport> processVCF(final String patient, final boolean isGermline, final File vcf,
            final BachelorEligibility eligibility)
    {
        final EligibilityReport.ReportType type =
                isGermline ? EligibilityReport.ReportType.GERMLINE_MUTATION : EligibilityReport.ReportType.SOMATIC_MUTATION;

        LOGGER.info("processing {} vcf: {}", type, vcf.getPath());

        try (final VCFFileReader reader = new VCFFileReader(vcf, true))
        {
            // assume that the first sample is the germline
            final String sample = reader.getFileHeader().getGenotypeSamples().get(0);
            return eligibility.processVCF(patient, sample, type, reader);
        }
        catch (final TribbleException e)
        {
            LOGGER.error("error with VCF file {}: {}", vcf.getPath(), e.getMessage());
            return Collections.emptyList();
        }
    }

    private static Collection<EligibilityReport> processPurpleCNV(final String patient, final File cnv,
            final BachelorEligibility eligibility) {
        LOGGER.info("processing cnv: {}", cnv.getPath());
        try
        {
            final List<GeneCopyNumber> copyNumbers = GeneCopyNumberFile.read(cnv);
            return eligibility.processCopyNumbers(patient, copyNumbers);
        }
        catch (final IOException e)
        {
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

    private void process(final BachelorEligibility eligibility, final RunDirectory run, final String sampleId)
    {
        final String patient = run.getPatientID();
        final boolean doGermline = run.germline() != null && mRunGermline;
        final boolean doSomatic = run.somatic() != null && mRunSomatic;
        final boolean doCopyNumber = run.copyNumber() != null && mRunCopyNumber;
        final boolean doStructuralVariants = run.structuralVariants() != null && mRunStructuralVariants;

        LOGGER.info("processing run for patient({}) from file({})", patient, run.prefix());

        final List<EligibilityReport> result = Lists.newArrayList();
        if (doGermline)
        {
            result.addAll(processVCF(patient, true, run.germline(), eligibility));
        }
        if (doSomatic)
        {
            result.addAll(processVCF(patient, false, run.somatic(), eligibility));
        }
        if (doCopyNumber)
        {
            result.addAll(processPurpleCNV(patient, run.copyNumber(), eligibility));
        }
        if (doStructuralVariants)
        {
            result.addAll(processSV(patient, run.structuralVariants(), eligibility));
        }

        if(result.isEmpty())
            return;

        createOutputFiles();

        try
        {
            for (final EligibilityReport r : result)
            {
                mMainDataWriter.write(String.format("%s,%s,%s,%s,%s,%s,%s,%d,%s,%s,%s,%s,%s,%s,%d,%s",
                        sampleId, r.source().toString(), r.program(), r.id(), r.genes(), r.transcriptId(),
                        r.chrom(), r.pos(), r.ref(), r.alts(), r.effects(), r.annotations(),
                        r.hgvsProtein(), r.isHomozygous(), r.phredScore(), r.hgvsCoding()));

                mMainDataWriter.newLine();

                mBedFileWriter.write(String.format("%s\t%s\t%d\t%d", sampleId, r.chrom(), r.pos() - 1, r.pos()));
                mBedFileWriter.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing output: {}", e.toString());
        }
    }

    private void createOutputFiles()
    {
        if(mMainDataWriter != null || mBedFileWriter != null)
            return;

        String outputDir = mOutputDir;

        if (!outputDir.endsWith("/"))
            outputDir += "/";

        String mainFileName = outputDir + "bachelor_output.csv";
        String bedFileName = outputDir + "bachelor_bed.csv";

        try
        {
            mMainDataWriter = Files.newBufferedWriter(Paths.get(mainFileName));
            mMainDataWriter.write(fileHeader());
            mMainDataWriter.newLine();

            mBedFileWriter = Files.newBufferedWriter(Paths.get(bedFileName));
        }
        catch(IOException e)
        {
            LOGGER.error("failed to create output files: {}", e.toString());
        }
    }

    private String fileHeader()
    {
        return String.join(",",
                Arrays.asList("SAMPLEID", "SOURCE", "PROGRAM", "ID", "GENE", "TRANSCRIPT_ID", "CHROM", "POS",
                        "REF", "ALTS", "EFFECTS",  "ANNOTATIONS", "HGVS_PROTEIN", "IS_HOMOZYGOUS", "PHRED_SCORE", "HGVS_CODING"));
    }

    public static void main(final String... args)
    {
        final Options options = createOptions();

        BachelorApplication bachelorApp = new BachelorApplication();

        try
        {
            final CommandLine cmd = createCommandLine(options, args);

            if (cmd.hasOption(LOG_DEBUG))
                Configurator.setRootLevel(Level.DEBUG);

            if(!bachelorApp.loadConfig(cmd))
            {
                System.exit(1);
                return;
            }

            if(!bachelorApp.run())
            {
                System.exit(1);
                return;
            }
        }
        catch (final ParseException e)
        {
            printHelpAndExit(options);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    private static void printHelpAndExit(final Options options)
    {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Bachelor", "Determines eligibility", options, "", true);
        System.exit(1);
    }

}

package com.hartwig.hmftools.bachelor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;

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

import htsjdk.tribble.TribbleException;
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
    private static final String SAMPLE = "sample";
    private static final String LOG_DEBUG = "log_debug";
    private static final String BATCH_MAX_DIR = "max_batch_dir"; // only for testingt

    private static final String BATCH_DIR_FILE = "batch_dir_file";

    private Map<String, Program> mProgramMap;
    private BufferedWriter mMainDataWriter;
    private BufferedWriter mBedFileWriter;

    // config
    private String mOutputDir;
    private String mBatchDirectory;
    private String mBatchDirectoryFile;
    private String mRunDirectoy;
    private boolean mIsBatchRun;
    private boolean mIsSingleRun;
    private String mSampleId;
    private int mMaxBatchDirectories;

    public BachelorApplication()
    {
        mProgramMap = null;
        mMainDataWriter = null;
        mBedFileWriter = null;

        mIsSingleRun = false;
        mIsBatchRun = false;
        mBatchDirectoryFile = "";
        mBatchDirectory = "";
        mOutputDir = "";
        mMaxBatchDirectories = 0;
    }

    @NotNull
    private static Options createOptions()
    {
        final Options options = new Options();

        options.addOption(CONFIG_DIRECTORY, true, "folder to find program XMLs");
        options.addOption(CONFIG_XML, true, "single config XML to run");
        options.addOption(OUTPUT_DIR, true, "output file");
        options.addOption(RUN_DIRECTORY, true, "the run directory to look for inputs");
        options.addOption(BATCH_DIRECTORY, true, "runs directory to batch process");
        options.addOption(BATCH_MAX_DIR, true, "Max batch directories to batch process");
        options.addOption(VALIDATE, false, "only validate the configs");
        options.addOption(GERMLINE, false, "process the germline file");
        options.addOption(BATCH_DIR_FILE, true, "Optional: list of directories to search");
        options.addOption(SOMATIC, false, "process the somatic file");
        options.addOption(SAMPLE, true, "sample id");
        options.addOption(LOG_DEBUG, false, "Sets log level to Debug, off by default");
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

        if(cmd.hasOption(BATCH_DIRECTORY))
        {
            mBatchDirectory = cmd.getOptionValue(BATCH_DIRECTORY);
            mBatchDirectoryFile = cmd.getOptionValue(BATCH_DIR_FILE);
            mMaxBatchDirectories = Integer.parseInt(cmd.getOptionValue(BATCH_MAX_DIR, "0"));
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
                List<RunDirectory> runDirectories = Lists.newArrayList();

                if(!mBatchDirectoryFile.isEmpty())
                {
                    LOGGER.debug("loading batch directories from file");

                    final List<String> batchDirectories = loadBatchDirectories(mBatchDirectory, mBatchDirectoryFile);
                    runDirectories = batchDirectories.stream()
                            .map(Paths::get)
                            .map(RunDirectory::new)
                            .collect(Collectors.toList());
                }
                else
                {
                    final Path root = Paths.get(mBatchDirectory);

                    try (final Stream<Path> stream = Files.walk(root, 1, FileVisitOption.FOLLOW_LINKS).parallel())
                    {
                        runDirectories = stream.filter(p -> p.toFile().isDirectory())
                                .filter(p -> !p.equals(root))
                                .map(RunDirectory::new)
                                .collect(Collectors.toList());
                    }
                    catch (Exception e)
                    {
                        LOGGER.error("failed walking batch directories: {}", e.toString());
                    }
                }

                LOGGER.info("found {} batch directories", runDirectories.size());

                // add the filtered and passed SV entries for each file
                for (int i = 0; i < runDirectories.size(); ++i)
                {
                    final RunDirectory runDir = runDirectories.get(i);

                    process(eligibility, runDir, "");

                    if(mMaxBatchDirectories > 0 && i >= mMaxBatchDirectories)
                        break;
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

    private static int BATCH_DIRECTORIES_CSV_FIELDS= 4;

    private static List<String> loadBatchDirectories(final String rootDir, final String filename)
    {
        List<String> batchDirectories = Lists.newArrayList();

        if (filename.isEmpty())
            return batchDirectories;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine(); // skip header

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                // CSV fields Month,Date,Time,Directory

                if (items.length < BATCH_DIRECTORIES_CSV_FIELDS)
                {
                    LOGGER.error("invalid CSV item length({})", items.length);
                    return batchDirectories;
                }

                final String directory = rootDir + items[3];

                batchDirectories.add(directory);
            }

            LOGGER.debug("loaded {} batch directories", batchDirectories.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read batch dirs input CSV file({})", filename);
            return batchDirectories;
        }

        return batchDirectories;
    }

    private static Collection<EligibilityReport> processVCF(final String patient, final String sampleId, final File vcf,
            final BachelorEligibility eligibility)
    {
        if(vcf == null)
            return Collections.emptyList();

        final EligibilityReport.ReportType type = EligibilityReport.ReportType.GERMLINE_MUTATION;

        LOGGER.debug("processing {} vcf: {}", type, vcf.getPath());

        try (final VCFFileReader reader = new VCFFileReader(vcf, true))
        {
            // assume that the first sample is the germline
            // final String sampleId = reader.getFileHeader().getGenotypeSamples().size() > 1 ? reader.getFileHeader().getGenotypeSamples().get(1) : reader.getFileHeader().getGenotypeSamples().get(0);
            return eligibility.processVCF(patient, sampleId, type, reader);
        }
        catch (final TribbleException e)
        {
            LOGGER.error("error with VCF file {}: {}", vcf.getPath(), e.getMessage());
            return Collections.emptyList();
        }
    }

    private void process(final BachelorEligibility eligibility, final RunDirectory runDir, String sampleId)
    {
        final String patient = runDir.getPatientID();

        if(sampleId.isEmpty())
        {
            try
            {
                final RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDir.sampleDir().toString());

                if(runContext != null)
                    sampleId = runContext.tumorSample();
            }
            catch (Exception e)
            {
                // skip using meta data
                sampleId = patient;
            }
        }

        LOGGER.info("processing run for patient({}) sampleId({}) directory({})", patient, sampleId, runDir.sampleDir());

        final List<EligibilityReport> result = Lists.newArrayList();
        result.addAll(processVCF(patient, sampleId, runDir.germline(), eligibility));

        if(result.isEmpty())
            return;

        createOutputFiles();

        try
        {
            for (final EligibilityReport r : result)
            {
                mMainDataWriter.write(String.format("%s,%s,%s,%s,%s,%s",
                        r.sampleId(), r.source().toString(), r.program(), r.id(), r.genes(), r.transcriptId()));

                mMainDataWriter.write(String.format(",%s,%d,%s,%s,%s,%s",
                        r.chrom(), r.pos(), r.ref(), r.alts(), r.effects(), r.annotations()));

                mMainDataWriter.write(String.format(",%s,%s,%d,%s,%s,%d,%d,%s",
                        r.hgvsProtein(), r.isHomozygous(), r.phredScore(), r.hgvsCoding(), r.matchType(), r.altCount(), r.readDepth(), r.condonInfo()));

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
            mMainDataWriter = Files.newBufferedWriter(Paths.get(mainFileName), StandardOpenOption.CREATE);

            mMainDataWriter.write("SAMPLEID,SOURCE,PROGRAM,ID,GENE,TRANSCRIPT_ID,CHROM,POS,REF,ALTS");
            mMainDataWriter.write(",EFFECTS,ANNOTATIONS,HGVS_PROTEIN,IS_HOMOZYGOUS,PHRED_SCORE,HGVS_CODING,MATCH_TYPE,ALT_COUNT,READ_DEPTH,CODON_INFO");

            mMainDataWriter.newLine();

            mBedFileWriter = Files.newBufferedWriter(Paths.get(bedFileName), StandardOpenOption.CREATE);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to create output files: {}", e.toString());
        }
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

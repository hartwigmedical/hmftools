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
    private static final String BATCH_MAX_DIR = "max_batch_dir"; // only for testing

    private static final String SAMPLE_LIST_FILE = "sample_list_file";

    private BachelorProgram mProgram;

    private Map<String, Program> mProgramMap;
    private BufferedWriter mMainDataWriter;
    private BufferedWriter mBedFileWriter;

    // config
    private String mOutputDir;
    private String mBatchDirectory;
    private String mRunDirectoy;
    private boolean mIsBatchRun;
    private boolean mIsSingleRun;
    private String mSampleId;
    private List<String> mLimitedSampleList;
    private int mMaxBatchDirectories;

    public BachelorApplication()
    {
        mProgramMap = null;
        mMainDataWriter = null;
        mBedFileWriter = null;

        mSampleId = "";
        mIsSingleRun = false;
        mIsBatchRun = false;
        mBatchDirectory = "";
        mOutputDir = "";
        mMaxBatchDirectories = 0;
        mLimitedSampleList = Lists.newArrayList();
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
        options.addOption(SAMPLE_LIST_FILE, true, "Optional: limiting list of sample IDs to process");
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
            LOGGER.error("no Programs loaded, exiting");
            return false;
        }

        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);

        if(cmd.hasOption(BATCH_DIRECTORY))
        {
            mBatchDirectory = cmd.getOptionValue(BATCH_DIRECTORY);
            mMaxBatchDirectories = Integer.parseInt(cmd.getOptionValue(BATCH_MAX_DIR, "0"));
            mIsBatchRun = true;

            if(cmd.hasOption(SAMPLE_LIST_FILE))
                loadSampleListFile(cmd.getOptionValue(SAMPLE_LIST_FILE));
        }
        else if(cmd.hasOption(RUN_DIRECTORY))
        {
            mSampleId = cmd.getOptionValue(SAMPLE);
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
        mProgram = new BachelorProgram();

        if(!mProgram.loadConfig(mProgramMap))
            return false;

        // final BachelorEligibility eligibility = BachelorEligibility.fromMap(mProgramMap);

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

                LOGGER.info("found {} batch directories", runDirectories.size());

                // add the filtered and passed SV entries for each file
                for (int i = 0; i < runDirectories.size(); ++i)
                {
                    final RunDirectory runDir = runDirectories.get(i);

                    processSampleDirectory(runDir, "");

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

                processSampleDirectory(new RunDirectory(path), mSampleId);
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

    private void loadSampleListFile(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine(); // skip header

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                final String sampleId = items[0];

                mLimitedSampleList.add(items[0]);
            }

            LOGGER.info("loaded {} specific sample ids", mLimitedSampleList.size());

        }
        catch (IOException exception)
        {
            LOGGER.error("Failed to read sample list input CSV file({}): {}", filename, exception.toString());
        }
    }

    private List<EligibilityReport> processVCF(final String sampleId, final File vcf)
    {
        if(vcf == null)
            return Lists.newArrayList();

        LOGGER.debug("processing vcf: {}", vcf.getPath());

        try (final VCFFileReader reader = new VCFFileReader(vcf, true))
        {
            return mProgram.processVCF(sampleId, reader);
        }
        catch (final TribbleException e)
        {
            LOGGER.error("error with VCF file {}: {}", vcf.getPath(), e.getMessage());
            return Lists.newArrayList();
        }
    }

    private void processSampleDirectory(final RunDirectory runDir, String sampleId)
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

        if(!mLimitedSampleList.isEmpty() && !mLimitedSampleList.contains(sampleId))
        {
            LOGGER.info("skipping sampleId({}) not in specified list", sampleId);
            return;
        }

        LOGGER.info("processing run for patient({}) sampleId({}) directory({})", patient, sampleId, runDir.sampleDir());

        final List<EligibilityReport> result = Lists.newArrayList();

        result.addAll(processVCF(sampleId, runDir.germline()));

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

                mMainDataWriter.write(String.format(",%s,%s,%d,%s,%s,%d,%d,%d,%d,%s",
                        r.hgvsProtein(), r.isHomozygous(), r.phredScore(), r.hgvsCoding(), r.matchType(),
                        r.germlineAltCount(), r.germlineReadDepth(), r.tumorAltCount(), r.tumorReadDepth(), r.condonInfo()));

                mMainDataWriter.newLine();

                if(!mIsBatchRun)
                {
                    mBedFileWriter.write(String.format("%s\t%s\t%d\t%d", sampleId, r.chrom(), r.pos() - 1, r.pos()));
                    mBedFileWriter.newLine();
                }
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
            mMainDataWriter.write(",EFFECTS,ANNOTATIONS,HGVS_PROTEIN,IS_HOMOZYGOUS,PHRED_SCORE,HGVS_CODING");
            mMainDataWriter.write(",MATCH_TYPE,GL_ALT_COUNT,GL_READ_DEPTH,TUMOR_ALT_COUNT,TUMOR_READ_DEPTH,CODON_INFO");

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

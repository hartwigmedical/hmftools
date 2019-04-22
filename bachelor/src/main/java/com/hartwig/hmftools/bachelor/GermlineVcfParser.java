package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.bachelor.types.ConfigSchema;
import com.hartwig.hmftools.bachelor.types.RunDirectory;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;

import nl.hartwigmedicalfoundation.bachelor.Program;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.xml.sax.SAXException;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.vcf.VCFFileReader;

public class GermlineVcfParser
{
    private Map<String, Program> mProgramMap;
    private BachelorProgram mProgram;
    ExternalDBFilters mFilterFileBuilder;

    // config
    private List<String> mSampleIds;
    private String mSampleDataDir;
    private String mOutputDir;
    private boolean mIsBatchMode;
    private boolean mSkipIndexFile;
    private String mSampleDataDirectory;
    private String mExternalFiltersFile;
    private int mMaxBatchDirectories;

    private BufferedWriter mVcfDataWriter;
    private BufferedWriter mBedFileWriter;

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);


    public GermlineVcfParser()
    {
        mProgramMap = null;
        mFilterFileBuilder = null;

        mSampleIds = Lists.newArrayList();
        mSampleDataDir = "";
        mOutputDir= "";
        mExternalFiltersFile = "";
        mSkipIndexFile = false;
        mVcfDataWriter = null;
        mBedFileWriter = null;
        mMaxBatchDirectories = 0;
    }

    private static final String CONFIG_XML = "configXml";
    private static final String SKIP_INDEX_FILE = "skip_index_file";
    private static final String EXTERNAL_FILTER_FILE = "ext_filter_file";
    private static final String CREATE_FILTER_FILE = "create_filter_file";
    private static final String BATCH_MAX_DIR = "max_batch_dir"; // only for testing

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(CONFIG_XML, true, "single config XML to run");
        options.addOption(BATCH_MAX_DIR, true, "Max batch directories to batch process");
        options.addOption(SKIP_INDEX_FILE, false, "Skip VCF index file");
        options.addOption(EXTERNAL_FILTER_FILE, true, "Optional: name of an external filter file");
        options.addOption(CREATE_FILTER_FILE, true, "Optional: create black and white list filter files");
    }

    public boolean initialise(final CommandLine cmd, final List<String> sampleIds, final String dataPath, final String outputDir)
    {
        try
        {
            if (cmd.hasOption(CONFIG_XML))
            {
                mProgramMap = loadXML(Paths.get(cmd.getOptionValue(CONFIG_XML)));
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

        if (mProgramMap.isEmpty())
        {
            LOGGER.error("no Programs loaded, exiting");
            return false;
        }

        if(cmd.hasOption(CREATE_FILTER_FILE))
        {
            LOGGER.info("building Clinvar filter files");
            final String filterInputFile = cmd.getOptionValue(CREATE_FILTER_FILE);
            mFilterFileBuilder = new ExternalDBFilters(filterInputFile);
            return true;
        }

        mExternalFiltersFile = cmd.getOptionValue(EXTERNAL_FILTER_FILE, "");
        mSkipIndexFile = cmd.hasOption(SKIP_INDEX_FILE);

        mIsBatchMode = (sampleIds.size() != 1);
        mMaxBatchDirectories = Integer.parseInt(cmd.getOptionValue(BATCH_MAX_DIR, "0"));

        mSampleIds.addAll(sampleIds);
        mSampleDataDir = dataPath;
        mOutputDir = outputDir;

        return true;
    }

    public boolean run()
    {
        if(mFilterFileBuilder != null)
        {
            final Program program = mProgramMap.values().iterator().next();

            mFilterFileBuilder.createFilterFile(mOutputDir, program);
            LOGGER.info("run complete");
            return true;
        }

        mProgram = new BachelorProgram();

        if(!mProgram.loadConfig(mProgramMap))
            return false;

        if(!mExternalFiltersFile.isEmpty())
        {
            mProgram.addExternalFilters(ExternalDBFilters.loadExternalFilters(mExternalFiltersFile));
        }

        if (mIsBatchMode)
        {
            LOGGER.info("beginning batch run");
        }
        else
        {
            LOGGER.info("beginning single sample run: {}", mSampleIds.get(0));
        }

        if (mIsBatchMode)
        {
            List<RunDirectory> runDirectories = Lists.newArrayList();

            final Path root = Paths.get(mSampleDataDirectory);

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
        else
        {
            final Path path = Paths.get(mSampleDataDir);

            if (!Files.exists(path))
            {
                LOGGER.error("-runDirectory path does not exist");
                return false;
            }

            processSampleDirectory(new RunDirectory(path), mSampleIds.get(0));
        }

        LOGGER.info("germline VCF parsing complete");

        closeBufferedWriter(mVcfDataWriter);
        closeBufferedWriter(mBedFileWriter);

        return true;
    }

    private void processSampleDirectory(final RunDirectory runDir, String sampleId)
    {
        final String patient = runDir.getPatientID();

        if(sampleId.isEmpty())
        {
            try
            {
                final RunContext runContext = ProductionRunContextFactory.fromRunDirectory(runDir.sampleDir().toString());
                sampleId = runContext.tumorSample();
            }
            catch (Exception e)
            {
                // Skip using meta data
                sampleId = patient;
            }
        }

        if(mSampleIds.size() > 1 && !mSampleIds.contains(sampleId))
        {
            LOGGER.info("skipping sampleId({}) not in specified list", sampleId);
            return;
        }

        LOGGER.info("Processing run for patient({}) sampleId({}) directory({})", patient, sampleId, runDir.sampleDir());

        final List<EligibilityReport> result = Lists.newArrayList();

        result.addAll(processVCF(sampleId, runDir.germline()));

        if(result.isEmpty())
            return;

        createOutputFiles();

        try
        {
            for (final EligibilityReport r : result)
            {
                mVcfDataWriter.write(String.format("%s,%s,%s,%s,%s",
                        r.sampleId(), r.program(), r.id(), r.genes(), r.transcriptId()));

                mVcfDataWriter.write(String.format(",%s,%d,%s,%s,%s,%s,%s",
                        r.chrom(), r.pos(), r.ref(), r.alts(), r.codingEffect(), r.effects(), r.annotations()));

                mVcfDataWriter.write(String.format(",%s,%s,%d,%s,%s,%s,%d,%d,%d,%d,%s",
                        r.hgvsProtein(), r.isHomozygous(), r.phredScore(), r.hgvsCoding(), r.matchType(),
                        r.hasDepthInfo(), r.germlineAltCount(), r.germlineReadDepth(), r.tumorAltCount(), r.tumorReadDepth(),
                        r.condonInfo()));

                mVcfDataWriter.write(String.format(",%s,%s,%s",
                        r.matchesClinvarFilter(), r.clinvarSignificance(), r.clinvarSigInfo()));

                mVcfDataWriter.newLine();

                if(!mIsBatchMode)
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

    private List<EligibilityReport> processVCF(final String sampleId, final File vcf)
    {
        if(vcf == null)
            return Lists.newArrayList();

        LOGGER.debug("Processing vcf: {}", vcf.getPath());

        try (final VCFFileReader reader = new VCFFileReader(vcf, !mSkipIndexFile))
        {
            return mProgram.processVcfFile(sampleId, reader, !mSkipIndexFile);
        }
        catch (final TribbleException e)
        {
            LOGGER.error("Error with VCF file {}: {}", vcf.getPath(), e.getMessage());
            return Lists.newArrayList();
        }
    }

    private void createOutputFiles()
    {
        if(mVcfDataWriter != null || mBedFileWriter != null)
            return;

        String outputDir = mOutputDir;

        if (!outputDir.endsWith("/"))
            outputDir += "/";

        String mainFileName = outputDir + "bachelor_output.csv";
        String bedFileName = outputDir + "bachelor_bed.csv";

        try
        {
            mVcfDataWriter = createBufferedWriter(mainFileName, false);

            mVcfDataWriter.write("SampleId,Program,Id,Gene,TranscriptId,Chromosome,Position,Ref,Alt");
            mVcfDataWriter.write(",CodingEffect,Effect,Annotations,HgvsProtein,IsHomozygous,PhredScore,HgvsCoding,");
            mVcfDataWriter.write(",MatchType,HasDepthInfo,GermlineAltCount,GermlineReadDepth,TumorAltCount,TumorReadDepth,CodonInfo");
            mVcfDataWriter.write(",ClinvarMatch,ClinvarSignificance,ClinvarSigInfo");


            mVcfDataWriter.newLine();

            mBedFileWriter = createBufferedWriter(bedFileName, false);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to create output files: {}", e.toString());
        }
    }

    private static Map<String, Program> loadXML(final Path path) throws IOException, SAXException
    {
        final ConfigSchema schema = ConfigSchema.make();

        final List<Program> programs = Files.walk(path)
                .filter(p -> p.toString().endsWith(".xml"))
                .map(schema::processXML)
                .filter(Objects::nonNull)
                .collect(Collectors.toList());

        final Map<String, Program> result = Maps.newHashMap();

        for (final Program p : programs)
        {
            if (result.containsKey(p.getName()))
            {
                LOGGER.error("duplicate Programs detected: {}", p.getName());
                System.exit(1);
            }
            else
            {
                result.put(p.getName(), p);
            }
        }

        return result;
    }

}

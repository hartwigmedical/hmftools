package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.BachelorApplication.DEFAULT_BACH_DIRECTORY;
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
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelor.types.RunDirectory;
import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
import com.hartwig.hmftools.common.context.RunContext;

import nl.hartwigmedicalfoundation.bachelor.Program;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.vcf.VCFFileReader;

public class GermlineVcfParser
{
    private Map<String, Program> mProgramConfigMap;
    private GermlineVariantFinder mProgram;

    // config
    private boolean mIsBatchMode;
    private String mBatchDataDir;
    private boolean mUsingBatchOutput;
    private boolean mSkipIndexFile;
    private String mExternalFiltersFile;

    private BufferedWriter mVcfDataWriter;
    private BufferedWriter mBedFileWriter;

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);


    public GermlineVcfParser()
    {
        mProgramConfigMap = null;
        mExternalFiltersFile = "";
        mBatchDataDir = "";
        mUsingBatchOutput = false;
        mSkipIndexFile = false;
        mVcfDataWriter = null;
        mBedFileWriter = null;
    }

    private static final String SKIP_INDEX_FILE = "skip_index_file";
    private static final String EXTERNAL_FILTER_FILE = "ext_filter_file";

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(SKIP_INDEX_FILE, false, "Skip VCF index file");
        options.addOption(EXTERNAL_FILTER_FILE, true, "Optional: name of an external filter file");
    }

    public boolean initialise(final CommandLine cmd, Map<String, Program> configMap, boolean isBatchMode, final String batchOutputDir)
    {
        mProgramConfigMap = configMap;

        if (mProgramConfigMap.isEmpty())
        {
            LOGGER.error("no Programs loaded, exiting");
            return false;
        }

        mExternalFiltersFile = cmd.getOptionValue(EXTERNAL_FILTER_FILE, "");
        mSkipIndexFile = cmd.hasOption(SKIP_INDEX_FILE);

        mBatchDataDir = batchOutputDir;
        mIsBatchMode = isBatchMode;
        mUsingBatchOutput = mIsBatchMode && !mBatchDataDir.isEmpty();

        mProgram = new GermlineVariantFinder();

        if(!mProgram.loadConfig(mProgramConfigMap))
            return false;

        if(!mExternalFiltersFile.isEmpty())
        {
            mProgram.addExternalFilters(ExternalDBFilters.loadExternalFilters(mExternalFiltersFile));
        }

        if(mUsingBatchOutput)
        {
            createOutputFiles(mBatchDataDir);
        }

        return true;
    }

    public List<BachelorGermlineVariant> getBachelorRecords() { return mProgram.getVariants(); }

    public void run(final RunDirectory runDir, String sampleId)
    {
        LOGGER.info("Processing run for sampleId({}) directory({})", sampleId, runDir.sampleDir());

        processVCF(sampleId, runDir.germline());

        if(mProgram.getVariants().isEmpty())
            return;

        if(!mUsingBatchOutput)
        {
            String sampleOutputDir = runDir.sampleDir().toString() + File.separator + DEFAULT_BACH_DIRECTORY + File.separator;
            createOutputFiles(sampleOutputDir);
        }

        try
        {
            for (final BachelorGermlineVariant variant : mProgram.getVariants())
            {
                mVcfDataWriter.write(variant.asCsv(true));
                mVcfDataWriter.newLine();

                if(mBedFileWriter != null)
                {
                    mBedFileWriter.write(String.format("%s\t%s\t%d\t%d",
                            sampleId, variant.Chromosome, variant.Position - 1, variant.Position));
                    mBedFileWriter.newLine();
                }
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing output: {}", e.toString());
        }

        if(!mUsingBatchOutput)
        {
            close();
        }
    }

    public void close()
    {
        if(mUsingBatchOutput)
        {
            closeBufferedWriter(mVcfDataWriter);
            mVcfDataWriter = null;
            closeBufferedWriter(mBedFileWriter);
            mBedFileWriter = null;
        }
    }

    private void processVCF(final String sampleId, final File vcf)
    {
        if(vcf == null)
            return;

        LOGGER.debug("Processing vcf: {}", vcf.getPath());

        try (final VCFFileReader reader = new VCFFileReader(vcf, !mSkipIndexFile))
        {
            mProgram.processVcfFile(sampleId, reader, !mSkipIndexFile);
        }
        catch (final TribbleException e)
        {
            LOGGER.error("Error with VCF file {}: {}", vcf.getPath(), e.getMessage());
        }
    }

    private void createOutputFiles(final String outputDir)
    {
        if(mVcfDataWriter != null || mBedFileWriter != null)
            return;

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

            if(!mIsBatchMode)
            {
                mBedFileWriter = createBufferedWriter(bedFileName, false);
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to create output files: {}", e.toString());
        }
    }

}

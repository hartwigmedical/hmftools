package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BATCH_FILE;
import static com.hartwig.hmftools.bachelor.types.BachelorConfig.INTERIM_FILENAME;
import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.bachelor.datamodel.Program;
import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.bachelor.types.RunDirectory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.vcf.VCFFileReader;

public class GermlineVcfParser
{
    private final GermlineVariantFinder mProgram;

    private final BachelorConfig mConfig;

    // config
    private final boolean mSkipIndexFile;
    private final String mExternalFiltersFile;

    private BufferedWriter mVcfDataWriter;
    private BufferedWriter mBedFileWriter;

    private static final Logger LOGGER = LogManager.getLogger(BachelorApplication.class);

    GermlineVcfParser(final BachelorConfig config, final CommandLine cmd)
    {
        mConfig = config;

        mVcfDataWriter = null;
        mBedFileWriter = null;

        mExternalFiltersFile = cmd.getOptionValue(EXTERNAL_FILTER_FILE, "");
        mSkipIndexFile = cmd.hasOption(SKIP_INDEX_FILE);

        mProgram = new GermlineVariantFinder();

        if(!mProgram.loadConfig(mConfig.ProgramConfigMap))
        {
            return;
        }

        if(!mExternalFiltersFile.isEmpty())
        {
            mProgram.addExternalFilters(ExternalDBFilters.loadExternalFilters(mExternalFiltersFile));
        }

        if(mConfig.IsBatchMode)
        {
            createOutputFiles(mConfig.OutputDir, BATCH_FILE);
        }
    }

    private static final String SKIP_INDEX_FILE = "skip_index_file";
    private static final String EXTERNAL_FILTER_FILE = "ext_filter_file";

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(SKIP_INDEX_FILE, false, "Skip VCF index file");
        options.addOption(EXTERNAL_FILTER_FILE, true, "Optional: name of an external filter file");
    }

    public List<BachelorGermlineVariant> getBachelorRecords() { return mProgram.getVariants(); }

    public void run(RunDirectory runDir, String sampleId, String singleSampleOutputDir)
    {
        if (mConfig.ProgramConfigMap.isEmpty())
        {
            LOGGER.error("No programs loaded, exiting");
            return;
        }

        LOGGER.info("Processing run for sampleId({}) directory({})", sampleId, runDir.sampleDir());

        processVCF(sampleId, runDir.germline());

        if(mProgram.getVariants().isEmpty())
        {
            LOGGER.debug("No valid variants found in " + runDir.germline().getPath());
            return;
        }

        if(!mConfig.IsBatchMode)
        {
            LOGGER.debug("Creating single sample output dir: " + singleSampleOutputDir);
            createOutputFiles(singleSampleOutputDir, sampleId);
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

        if(!mConfig.IsBatchMode)
        {
            close();
        }
    }

    public void close()
    {
        closeBufferedWriter(mVcfDataWriter);
        mVcfDataWriter = null;
        closeBufferedWriter(mBedFileWriter);
        mBedFileWriter = null;
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

    private void createOutputFiles(final String outputDir, final String filePrefix)
    {
        if(mVcfDataWriter != null || mBedFileWriter != null)
            return;

        final String mainFileName = outputDir + filePrefix + INTERIM_FILENAME;

        try
        {
            mVcfDataWriter = createBufferedWriter(mainFileName, false);

            mVcfDataWriter.write("SampleId,Program,Id,Gene,TranscriptId,Chromosome,Position,Ref,Alt");
            mVcfDataWriter.write(",CodingEffect,Effect,Annotations,HgvsProtein,IsHomozygous,PhredScore,HgvsCoding,");
            mVcfDataWriter.write(",MatchType,HasDepthInfo,GermlineAltCount,GermlineReadDepth,TumorAltCount,TumorReadDepth,CodonInfo");
            mVcfDataWriter.write(",ClinvarMatch,ClinvarSignificance,ClinvarSigInfo");


            mVcfDataWriter.newLine();

            if(!mConfig.IsBatchMode)
            {
                String bedFileName = outputDir + "bachelor_bed.csv";
                mBedFileWriter = createBufferedWriter(bedFileName, false);
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to create output files: {}", e.toString());
        }
    }
}

package com.hartwig.hmftools.bachelor;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;

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

    private static final Logger LOGGER = LogManager.getLogger(GermlineVcfParser.class);

    GermlineVcfParser(final BachelorConfig config, final CommandLine cmd)
    {
        mConfig = config;

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
    }

    private static final String SKIP_INDEX_FILE = "skip_index_file";
    private static final String EXTERNAL_FILTER_FILE = "ext_filter_file";

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(SKIP_INDEX_FILE, false, "Skip VCF index file");
        options.addOption(EXTERNAL_FILTER_FILE, true, "Optional: name of an external filter file");
    }

    public List<BachelorGermlineVariant> getBachelorRecords() { return mProgram.getVariants(); }

    public void run(final String vcfFile, String sampleId, String singleSampleOutputDir)
    {
        if (mConfig.ProgramConfigMap.isEmpty())
        {
            LOGGER.error("No programs loaded, exiting");
            return;
        }

        if(!Files.exists(Paths.get(vcfFile)))
        {
            LOGGER.info("sampleId({}) germline VCF({}) not found", sampleId, vcfFile);
            return;
        }

        LOGGER.info("sampleId({}) reading germline VCF({})", sampleId, vcfFile);

        final File germlineVcf = new File(vcfFile);

        processVCF(sampleId, germlineVcf);

        if(mProgram.getVariants().isEmpty())
        {
            LOGGER.debug("no valid variants found");
            return;
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
}

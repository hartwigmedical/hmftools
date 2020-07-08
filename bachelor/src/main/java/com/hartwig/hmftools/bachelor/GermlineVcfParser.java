package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.bachelor.types.BachelorConfig;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.vcf.VCFFileReader;

public class GermlineVcfParser
{
    private final BachelorConfig mConfig;
    private final GermlineVariantFinder mProgram;
    private final boolean mSkipIndexFile;

    GermlineVcfParser(final BachelorConfig config, final CommandLine cmd)
    {
        mConfig = config;

        mProgram = new GermlineVariantFinder(mConfig.IncludeVcfFiltered);

        mSkipIndexFile = cmd.hasOption(SKIP_INDEX_FILE);

        if(!mProgram.loadConfig(mConfig.ProgramConfigMap))
        {
            return;
        }

        final String externalFiltersFile = cmd.getOptionValue(EXTERNAL_FILTER_FILE, "");
        if(!externalFiltersFile.isEmpty())
        {
            mProgram.addExternalFilters(ExternalDBFilters.loadExternalFilters(externalFiltersFile));
        }
    }

    private static final String SKIP_INDEX_FILE = "skip_index_file";
    private static final String EXTERNAL_FILTER_FILE = "ext_filter_file";

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(SKIP_INDEX_FILE, false, "Skip VCF index file");
        options.addOption(EXTERNAL_FILTER_FILE, true, "Optional: name of an external filter file");
    }

    List<BachelorGermlineVariant> getBachelorRecords() { return mProgram.getVariants(); }

    boolean run(final String vcfFile, String sampleId)
    {
        if (mConfig.ProgramConfigMap.isEmpty())
        {
            BACH_LOGGER.error("No programs loaded, exiting");
            return false;
        }

        if(!Files.exists(Paths.get(vcfFile)))
        {
            BACH_LOGGER.info("SampleId({}) germline VCF({}) not found", sampleId, vcfFile);
            return false;
        }

        BACH_LOGGER.info("SampleId({}) reading germline VCF({})", sampleId, vcfFile);

        final File germlineVcf = new File(vcfFile);

        if(!processVCF(sampleId, germlineVcf))
            return false;

        if(mProgram.getVariants().isEmpty())
        {
            BACH_LOGGER.debug("No valid variants found");
        }

        return true;
    }

    private boolean processVCF(final String sampleId, final File vcf)
    {
        if(vcf == null)
            return false;

        BACH_LOGGER.debug("Processing vcf: {}", vcf.getPath());

        try (final VCFFileReader reader = new VCFFileReader(vcf, !mSkipIndexFile))
        {
            mProgram.processVcfFile(sampleId, reader, !mSkipIndexFile);
        }
        catch (final TribbleException e)
        {
            BACH_LOGGER.error("Error with VCF file {}: {}", vcf.getPath(), e.getMessage());
            return false;
        }

        return true;
    }
}

package com.hartwig.hmftools.purple.config;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

import com.hartwig.hmftools.common.variant.GenotypeIds;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.variant.vcf.VCFFileReader;

public class SampleDataFiles
{
    public final String SampleDataDir;
    public final String SomaticSvVcfFile;
    public final String GermlineSvVcfFile;
    public final String RecoveredSvVcfFile;
    public final String SomaticVcfFile;
    public final String GermlineVcfFile;
    public final String AmberDirectory;
    public final String CobaltDirectory;

    public static final String SAMPLE_DIR = "sample_dir";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";
    private static String STRUCTURAL_VARIANTS = "structural_vcf";
    private static String SOMATIC_SV_VCF = "somatic_sv_vcf";
    private static String GERMLINE_SV_VCF = "germline_sv_vcf";
    private static String SV_RECOVERY_VCF = "sv_recovery_vcf";
    public static String GERMLINE_VARIANTS = "germline_vcf";
    private static String SOMATIC_VARIANTS = "somatic_vcf";

    static void addOptions(final Options options)
    {
        options.addOption(SAMPLE_DIR, true,"Path to the sample's directory where expect to find Cobalt, Amber, Gripss etc directories");

        options.addOption(COBALT,true,
                "Path to COBALT output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/cobalt");

        options.addOption(AMBER, true,
                "Path to AMBER output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/amber");

        options.addOption(STRUCTURAL_VARIANTS, true, "(Deprecated for 'somatic_sv_vcf', location of somatic SV VCF");
        options.addOption(SOMATIC_SV_VCF, true, "Optional location of somatic SV VCF");
        options.addOption(GERMLINE_SV_VCF, true, "Optional location of germline SV VCF to annotate");
        options.addOption(SV_RECOVERY_VCF, true, "Optional location of failing structural variants that may be recovered");
        options.addOption(GERMLINE_VARIANTS, true, "Optional location of germline variants to enrich and process in driver catalog.");
        options.addOption(SOMATIC_VARIANTS, true, "Optional location of somatic variant vcf to assist fitting in highly-diploid samples.");
    }

    public SampleDataFiles(final CommandLine cmd, final String sampleId, final String referenceId) throws ParseException
    {
        SampleDataDir = cmd.hasOption(SAMPLE_DIR) ? checkAddDirSeparator(cmd.getOptionValue(SAMPLE_DIR)) : null;

        if(cmd.hasOption(AMBER))
            AmberDirectory = cmd.getOptionValue(AMBER);
        else if(SampleDataDir != null)
            AmberDirectory = SampleDataDir + "amber/";
        else
            throw new ParseException("missing amber or sample_data_dir config");

        if(cmd.hasOption(COBALT))
            CobaltDirectory = cmd.getOptionValue(COBALT);
        else if(SampleDataDir != null)
            CobaltDirectory = SampleDataDir + "cobalt/";
        else
            throw new ParseException("missing cobalt or sample_data_dir config");

        if(cmd.hasOption(STRUCTURAL_VARIANTS))
            SomaticSvVcfFile = getFilename(cmd, STRUCTURAL_VARIANTS, SampleDataDir, sampleId, ".gripss.filtered.somatic.vcf.gz");
        else
            SomaticSvVcfFile = getFilename(cmd, SOMATIC_SV_VCF, SampleDataDir, sampleId, ".gripss.filtered.somatic.vcf.gz");

        GermlineSvVcfFile = getFilename(cmd, GERMLINE_SV_VCF, SampleDataDir, sampleId, ".gripss.filtered.germline.vcf.gz");

        RecoveredSvVcfFile = getFilename(cmd, SV_RECOVERY_VCF, SampleDataDir, sampleId, ".gripss.somatic.vcf.gz");
        SomaticVcfFile = getFilename(cmd, SOMATIC_VARIANTS, SampleDataDir, sampleId, ".sage.somatic.filtered.pave.vcf.gz");
        GermlineVcfFile = getFilename(cmd, GERMLINE_VARIANTS, SampleDataDir, sampleId, ".sage.germline.filtered.pave.vcf.gz");
    }

    public boolean hasValidSampleNames(final PurpleConfig config)
    {
        return hasValidVcfSampleNames(SomaticSvVcfFile, config, false)
            && hasValidVcfSampleNames(GermlineSvVcfFile, config, true)
            && hasValidVcfSampleNames(SomaticVcfFile, config, false)
            && hasValidVcfSampleNames(GermlineVcfFile, config, true);
    }

    private boolean hasValidVcfSampleNames(final String vcfFile, final PurpleConfig config, boolean isGermline)
    {
        if(vcfFile.isEmpty())
            return true;

        boolean referenceFirst = !isGermline;

        VCFFileReader vcfReader = new VCFFileReader(new File(vcfFile), false);

        boolean validVcfNames = GenotypeIds.hasValidSampleIds(
                vcfReader.getFileHeader(), config.ReferenceId, config.TumorId,  referenceFirst, false);

        if(!validVcfNames)
        {
            PPL_LOGGER.error("vcf({}) has invalid sample names: {}", vcfFile, vcfReader.getFileHeader().getGenotypeSamples());
            return false;
        }

        return true;
    }

    private String getFilename(
            final CommandLine cmd, final String config, final String sampleDataDir,
            final String sampleId, final String fileSuffix) throws ParseException
    {
        if(cmd.hasOption(config))
        {
            final String filename = cmd.getOptionValue(config);

            if(Files.exists(Paths.get(filename)))
                return filename;
            else
                throw new ParseException(String.format("missing file: %s", filename));
        }

        if(sampleDataDir == null)
            return "";

        final String filename = sampleDataDir + sampleId + fileSuffix;

        return Files.exists(Paths.get(filename)) ? filename : "";
    }
}

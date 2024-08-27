package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.AMBER_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.COBALT_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.ESVEE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.GRIPSS_GERMLINE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.GRIPSS_SOMATIC_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PAVE_GERMLINE_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PAVE_SOMATIC_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.VcfFileReader;

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
    private static String SOMATIC_SV_VCF = "somatic_sv_vcf";
    private static String GERMLINE_SV_VCF = "germline_sv_vcf";
    private static String SV_RECOVERY_VCF = "sv_recovery_vcf";
    public static String GERMLINE_VARIANTS = "germline_vcf";
    private static String SOMATIC_VARIANTS = "somatic_vcf";

    static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(SAMPLE_DIR, false,"Path to the sample's directory where expect to find Cobalt, Amber, Gripss etc directories");

        configBuilder.addPath(COBALT,false,
                "Path to COBALT output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/cobalt");

        configBuilder.addPath(AMBER, false,
                "Path to AMBER output directory. Required if <run_dir> not set, otherwise defaults to <run_dir>/amber");

        configBuilder.addPath(SOMATIC_SV_VCF, false, "Optional location of somatic SV VCF");
        configBuilder.addPath(GERMLINE_SV_VCF, false, "Optional location of germline SV VCF to annotate");
        configBuilder.addPath(SV_RECOVERY_VCF, false, "Optional location of failing structural variants that may be recovered");
        configBuilder.addPath(GERMLINE_VARIANTS, false, "Optional location of germline variants to enrich and process in driver catalog.");
        configBuilder.addPath(SOMATIC_VARIANTS, false, "Optional location of somatic variant vcf to assist fitting in highly-diploid samples.");
    }

    public SampleDataFiles(final ConfigBuilder configBuilder, final String sampleId)
    {
        SampleDataDir = configBuilder.hasValue(SAMPLE_DIR) ? checkAddDirSeparator(configBuilder.getValue(SAMPLE_DIR)) : null;

        if(configBuilder.hasValue(AMBER))
            AmberDirectory = checkAddDirSeparator(configBuilder.getValue(AMBER));
        else if(SampleDataDir != null)
            AmberDirectory = SampleDataDir + AMBER_DIR + File.separator;
        else
            AmberDirectory = null;

        if(configBuilder.hasValue(COBALT))
            CobaltDirectory = checkAddDirSeparator(configBuilder.getValue(COBALT));
        else if(SampleDataDir != null)
            CobaltDirectory = SampleDataDir + COBALT_DIR + File.separator;
        else
            CobaltDirectory = null;

        String somaticGripssVcfFile = getFilename(configBuilder, SOMATIC_SV_VCF, GRIPSS_SOMATIC_DIR, sampleId, ".gripss.filtered.somatic.vcf.gz");
        String germlineGripssVcfFile = getFilename(configBuilder, GERMLINE_SV_VCF, GRIPSS_GERMLINE_DIR, sampleId, ".gripss.filtered.germline.vcf.gz");

        String somaticEsveeVcfFile = getFilename(configBuilder, SOMATIC_SV_VCF, ESVEE_DIR, sampleId, ".esvee.somatic.vcf.gz");
        String germlineEsveeVcfFile = getFilename(configBuilder, GERMLINE_SV_VCF, ESVEE_DIR, sampleId, ".esvee.germline.vcf.gz");

        if(Files.exists(Paths.get(somaticGripssVcfFile)) | Files.exists(Paths.get(germlineGripssVcfFile)))
        {
            SomaticSvVcfFile = somaticGripssVcfFile;
            GermlineSvVcfFile = germlineGripssVcfFile;
        }
        else
        {
            SomaticSvVcfFile = somaticEsveeVcfFile;
            GermlineSvVcfFile = germlineEsveeVcfFile;
        }

        RecoveredSvVcfFile = getFilename(configBuilder, SV_RECOVERY_VCF, GRIPSS_SOMATIC_DIR, sampleId, ".gripss.somatic.vcf.gz");

        SomaticVcfFile = getFilename(configBuilder, SOMATIC_VARIANTS, PAVE_SOMATIC_DIR, sampleId, ".pave.somatic.vcf.gz");
        GermlineVcfFile = getFilename(configBuilder, GERMLINE_VARIANTS, PAVE_GERMLINE_DIR, sampleId, ".pave.germline.vcf.gz");
    }

    public boolean usesGridssSVs() { return SomaticSvVcfFile.contains("gridss") || GermlineSvVcfFile.contains("gridss"); }

    public boolean hasValidSampleNames(final PurpleConfig config)
    {
        boolean germlineValid = !config.runGermline()
                || (hasValidVcfSampleNames(GermlineSvVcfFile, config) && hasValidVcfSampleNames(GermlineVcfFile, config));

        boolean tumorValid = !config.runGermline()
                || (hasValidVcfSampleNames(SomaticSvVcfFile, config) && hasValidVcfSampleNames(SomaticVcfFile, config));

        return tumorValid && germlineValid;
    }

    private boolean hasValidVcfSampleNames(final String vcfFile, final PurpleConfig config)
    {
        if(vcfFile.isEmpty())
            return true;

        VcfFileReader vcfReader = new VcfFileReader(vcfFile);

        List<String> vcfSampleNames = vcfReader.vcfHeader().getGenotypeSamples();

        if(config.runTumor() && !vcfSampleNames.contains(config.TumorId))
        {
            PPL_LOGGER.error("vcf({}) missing tumorId({}) genotype", vcfFile, config.TumorId);
            return false;
        }

        if(config.runGermline() && !vcfSampleNames.contains(config.ReferenceId))
        {
            PPL_LOGGER.error("vcf({}) missing referenceId({}) genotype", vcfFile, config.ReferenceId);
            return false;
        }

        return true;
    }

    private String getFilename(
            final ConfigBuilder configBuilder, final String config, final String toolDir, final String sampleId, final String fileSuffix)
    {
        if(configBuilder.hasValue(config))
        {
            final String filename = configBuilder.getValue(config);

            if(Files.exists(Paths.get(filename)))
                return filename;

            PPL_LOGGER.error("missing file: {}", filename);
            return null;
        }

        if(SampleDataDir == null)
            return "";

        String filename = SampleDataDir + toolDir + File.separator + sampleId + fileSuffix;

        if(Files.exists(Paths.get(filename)))
            return filename;

        filename = SampleDataDir + sampleId + fileSuffix;

        return Files.exists(Paths.get(filename)) ? filename : "";
    }
}

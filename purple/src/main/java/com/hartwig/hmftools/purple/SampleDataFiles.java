package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.AMBER_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.COBALT_DIR;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.ESVEE_DIR;
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
    public final String SomaticVcfFile;
    public final String GermlineVcfFile;
    public final String AmberDirectory;
    public final String CobaltDirectory;

    public static final String SAMPLE_DIR = "sample_dir";
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";
    private static String SOMATIC_SV_VCF = "somatic_sv_vcf";
    private static String GERMLINE_SV_VCF = "germline_sv_vcf";
    public static String GERMLINE_VARIANTS = "germline_vcf";
    private static String SOMATIC_VARIANTS = "somatic_vcf";

    static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(SAMPLE_DIR, false,
                "Sample data directory containing Cobalt, Amber, Pave and Esvee files");

        configBuilder.addPath(COBALT,false,
                "Cobalt directory - required if <run_dir> not set, otherwise defaults to <run_dir>/cobalt");

        configBuilder.addPath(AMBER, false,
                "Amber directory - required if <run_dir> not set, otherwise defaults to <run_dir>/amber");

        configBuilder.addPath(SOMATIC_SV_VCF, false, "Somatic SV VCF");
        configBuilder.addPath(GERMLINE_SV_VCF, false, "Germline SV VCF to annotate");
        configBuilder.addPath(GERMLINE_VARIANTS, false, "Germline variant VCF");
        configBuilder.addPath(SOMATIC_VARIANTS, false, "Somatic variant VCF");
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

        SomaticSvVcfFile = getFilename(configBuilder, SOMATIC_SV_VCF, ESVEE_DIR, sampleId, ".esvee.somatic.vcf.gz");;
        GermlineSvVcfFile = getFilename(configBuilder, GERMLINE_SV_VCF, ESVEE_DIR, sampleId, ".esvee.germline.vcf.gz");;
        SomaticVcfFile = getFilename(configBuilder, SOMATIC_VARIANTS, PAVE_SOMATIC_DIR, sampleId, ".pave.somatic.vcf.gz");
        GermlineVcfFile = getFilename(configBuilder, GERMLINE_VARIANTS, PAVE_GERMLINE_DIR, sampleId, ".pave.germline.vcf.gz");
    }

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

package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_CFG;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_FILE_CFG;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectoriesFile.COBALT_DIR;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.AMBER_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.AMBER_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ESVEE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ESVEE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PAVE_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PAVE_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PAVE_SOMATIC_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PAVE_SOMATIC_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_TUMOR_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REDUX_TUMOR_DIR_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
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
    public final String ReduxTumorDirectory;

    public static final String SAMPLE_DIR = "sample_dir";
    public static final String SOMATIC_SV_VCF = "somatic_sv_vcf";
    public static final String GERMLINE_SV_VCF = "germline_sv_vcf";
    public static final String GERMLINE_VARIANTS = "germline_vcf";
    public static final String SOMATIC_VARIANTS = "somatic_vcf";

    @Deprecated
    private static final String AMBER = "amber";
    private static final String COBALT = "cobalt";

    static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(SAMPLE_DIR, false,
                "Sample data directory containing Cobalt, Amber, Pave and Esvee files");

        configBuilder.addPath(COBALT_DIR_CFG, false, COBALT_DIR_DESC);
        configBuilder.addPath(COBALT, false, COBALT_DIR_DESC + " deprecated, use: " + COBALT_DIR_CFG);

        configBuilder.addPath(AMBER, false, AMBER_DIR_DESC + " deprecated, use: " + AMBER_DIR_CFG);
        configBuilder.addPath(AMBER_DIR_CFG, false, AMBER_DIR_DESC);

        configBuilder.addPath(ESVEE_DIR_CFG, false, ESVEE_DIR_DESC);
        configBuilder.addPath(PAVE_SOMATIC_DIR_CFG, false, PAVE_SOMATIC_DIR_DESC);
        configBuilder.addPath(PAVE_GERMLINE_DIR_CFG, false, PAVE_GERMLINE_DIR_DESC);

        configBuilder.addPath(SOMATIC_SV_VCF, false, "Somatic SV VCF");
        configBuilder.addPath(GERMLINE_SV_VCF, false, "Germline SV VCF to annotate");
        configBuilder.addPath(GERMLINE_VARIANTS, false, "Germline variant VCF");
        configBuilder.addPath(SOMATIC_VARIANTS, false, "Somatic variant VCF");
        configBuilder.addPath(REDUX_TUMOR_DIR_CFG, false, REDUX_TUMOR_DIR_DESC);

        PipelineToolDirectories.addPipelineFormatOptions(configBuilder);
    }

    public SampleDataFiles(final ConfigBuilder configBuilder, final String sampleId)
    {
        SampleDataDir = configBuilder.hasValue(SAMPLE_DIR) ? checkAddDirSeparator(configBuilder.getValue(SAMPLE_DIR)) : null;

        PipelineToolDirectories pipelineToolDirectories = PipelineToolDirectories.resolveToolDirectories(
                configBuilder, PIPELINE_FORMAT_CFG, PIPELINE_FORMAT_FILE_CFG);

        if(configBuilder.hasValue(AMBER_DIR_CFG))
        {
            AmberDirectory = checkAddDirSeparator(configBuilder.getValue(AMBER_DIR_CFG));
        }
        else if(configBuilder.hasValue(AMBER))
        {
            AmberDirectory = checkAddDirSeparator(configBuilder.getValue(AMBER));
        }
        else if(SampleDataDir != null)
        {
            AmberDirectory = SampleDataDir + pipelineToolDirectories.amberDir() + File.separator;
        }
        else
        {
            AmberDirectory = null;
        }

        if(configBuilder.hasValue(COBALT_DIR_CFG))
        {
            CobaltDirectory = checkAddDirSeparator(configBuilder.getValue(COBALT_DIR_CFG));
        }
        else if(configBuilder.hasValue(COBALT))
        {
            CobaltDirectory = checkAddDirSeparator(configBuilder.getValue(COBALT));
        }
        else if(SampleDataDir != null)
        {
            CobaltDirectory = SampleDataDir + pipelineToolDirectories.cobaltDir() + File.separator;
        }
        else
        {
            CobaltDirectory = null;
        }

        SomaticSvVcfFile = getFilename(configBuilder, SOMATIC_SV_VCF, ESVEE_DIR_CFG,
                pipelineToolDirectories.esveeDir(), sampleId, ".esvee.somatic.vcf.gz");

        GermlineSvVcfFile = getFilename(configBuilder, GERMLINE_SV_VCF, ESVEE_DIR_CFG,
                pipelineToolDirectories.esveeDir(), sampleId, ".esvee.germline.vcf.gz");

        SomaticVcfFile = getFilename(configBuilder, SOMATIC_VARIANTS, PAVE_SOMATIC_DIR_CFG,
                pipelineToolDirectories.paveSomaticDir(), sampleId, ".pave.somatic.vcf.gz");

        GermlineVcfFile = getFilename(configBuilder, GERMLINE_VARIANTS, PAVE_GERMLINE_DIR_CFG,
                pipelineToolDirectories.paveGermlineDir(), sampleId, ".pave.germline.vcf.gz");

        ReduxTumorDirectory = configBuilder.getValue(REDUX_TUMOR_DIR_CFG);
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
        {
            return true;
        }

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
            final ConfigBuilder configBuilder, final String filePathConfig, final String toolDirConfig,
            final String pipelineToolDir, final String sampleId, final String fileSuffix)
    {
        if(configBuilder.hasValue(filePathConfig))
        {
            final String filename = configBuilder.getValue(filePathConfig);

            if(Files.exists(Paths.get(filename)))
            {
                return filename;
            }

            PPL_LOGGER.error("missing file: {}", filename);
            return null;
        }

        final String baseDir;
        if(configBuilder.hasValue(toolDirConfig))
        {
            baseDir = checkAddDirSeparator(configBuilder.getValue(toolDirConfig));
        }
        else if(SampleDataDir != null)
        {
            baseDir = checkAddDirSeparator(SampleDataDir) + checkAddDirSeparator(pipelineToolDir);
        }
        else
        {
            return "";
        }

        String filename = baseDir + sampleId + fileSuffix;

        if(Files.exists(Paths.get(filename)))
            return filename;

        // handle earlier pave formats, for conveniences only
        if(pipelineToolDir.contains(GERMLINE_SUB_DIR) && Files.exists(Paths.get(filename.replaceFirst(GERMLINE_SUB_DIR, ""))))
            return filename.replaceFirst(GERMLINE_SUB_DIR, "");

        if(pipelineToolDir.contains(SOMATIC_SUB_DIR) && Files.exists(Paths.get(filename.replaceFirst(SOMATIC_SUB_DIR, ""))))
            return filename.replaceFirst(SOMATIC_SUB_DIR, "");

        filename = SampleDataDir + sampleId + fileSuffix;

        return Files.exists(Paths.get(filename)) ? filename : "";
    }

    private static String GERMLINE_SUB_DIR = "/germline";
    private static String SOMATIC_SUB_DIR = "/somatic";
}

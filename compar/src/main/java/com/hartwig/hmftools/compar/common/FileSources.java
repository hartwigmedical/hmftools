package com.hartwig.hmftools.compar.common;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_CFG;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_DESC;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_FILE_CFG;
import static com.hartwig.hmftools.common.pipeline.PipelineToolDirectories.PIPELINE_FORMAT_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CHORD_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CIDER_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CIDER_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.COBALT_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.CUPPA_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LILAC_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_GERMLINE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SIGS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TEAL_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TEAL_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.V_CHORD_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.V_CHORD_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FileSources
{
    public final SourceType Source;
    public final String Linx;
    public final String LinxGermline;
    public final String Cobalt;
    public final String Purple;
    public final String Cuppa;
    public final String Lilac;
    public final String Chord;
    public final String Peach;
    public final String Virus;
    public final String SomaticVcf;
    public final String SomaticUnfilteredVcf;
    public final String TumorFlagstat;
    public final String GermlineFlagstat;
    public final String TumorBamMetrics;
    public final String GermlineBamMetrics;
    public final String SnpGenotype;
    public final String Cider;
    public final String Teal;
    public final String VChord;
    public final String Sigs;
    public final String Isofox;

    private static final String SAMPLE_DIR = "sample_dir";
    private static final String SOMATIC_VCF = "somatic_vcf";
    private static final String SOMATIC_UNFILTERED_VCF = "somatic_unfiltered_vcf";
    private static final String TUMOR_FLAGSTAT = "tumor_flagstat_dir";
    private static final String GERMLINE_FLAGSTAT = "germline_flagstat_dir";
    private static final String TUMOR_BAM_METRICS = "tumor_bam_metrics_dir";
    private static final String GERMLINE_BAM_METRICS = "germline_bam_metrics_dir";
    private static final String SNP_GENOTYPE = "snp_genotype_dir";

    public FileSources(
            final SourceType source, final String linx, final String cobalt, final String purple, final String linxGermline, final String cuppa,
            final String lilac, final String chord, final String peach, final String virus, final String somaticVcf,
            final String somaticUnfilteredVcf, final String tumorFlagstat, final String germlineFlagstat, final String tumorBamMetrics,
            final String germlineBamMetrics, final String snpGenotype, final String cider, final String teal, final String vChord,
            final String sigs, final String isofox)
    {
        Source = source;
        Linx = linx;
        LinxGermline = linxGermline;
        Cobalt = cobalt;
        Purple = purple;
        Cuppa = cuppa;
        Lilac = lilac;
        Chord = chord;
        Peach = peach;
        Virus = virus;
        SomaticVcf = somaticVcf;
        SomaticUnfilteredVcf = somaticUnfilteredVcf;
        TumorFlagstat = tumorFlagstat;
        GermlineFlagstat = germlineFlagstat;
        TumorBamMetrics = tumorBamMetrics;
        GermlineBamMetrics = germlineBamMetrics;
        SnpGenotype = snpGenotype;
        Cider = cider;
        Teal = teal;
        VChord = vChord;
        Sigs = sigs;
        Isofox = isofox;
    }

    public static FileSources sampleInstance(final FileSources fileSources, final String sampleId, final String referenceId)
    {
        return new FileSources(
                fileSources.Source,
                convertWildcardSamplePath(fileSources.Linx, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Cobalt, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Purple, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.LinxGermline, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Cuppa, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Lilac, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Chord, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Peach, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Virus, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.SomaticVcf, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.SomaticUnfilteredVcf, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.TumorFlagstat, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.GermlineFlagstat, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.TumorBamMetrics, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.GermlineBamMetrics, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.SnpGenotype, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Cider, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Teal, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.VChord, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Sigs, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.Isofox, sampleId, referenceId));
    }

    private static void addPathConfig(
            final ConfigBuilder configBuilder, final String toolDir, final String toolDesc, final SourceType sourceType)
    {
        configBuilder.addPrefixedPath(
                formSourceConfig(toolDir, sourceType), false, formSourceDescription(toolDesc, sourceType),
                formSourceConfig(SAMPLE_DIR, sourceType));
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        for(SourceType sourceType : SourceType.values())
        {
            configBuilder.addPath(
                    formSourceConfig(SAMPLE_DIR, sourceType), false,
                    formSourceDescription("Sample data root directory", sourceType));

            addPathConfig(configBuilder, LINX_DIR_CFG, LINX_DIR_DESC, sourceType);
            addPathConfig(configBuilder, LINX_GERMLINE_DIR_CFG, LINX_GERMLINE_DIR_DESC, sourceType);
            addPathConfig(configBuilder, COBALT_DIR_CFG, COBALT_DIR_DESC, sourceType);
            addPathConfig(configBuilder, PURPLE_DIR_CFG, PURPLE_DIR_DESC, sourceType);
            addPathConfig(configBuilder, LILAC_DIR_CFG, LILAC_DIR_DESC, sourceType);
            addPathConfig(configBuilder, CHORD_DIR_CFG, CHORD_DIR_DESC, sourceType);
            addPathConfig(configBuilder, CUPPA_DIR_CFG, CUPPA_DIR_DESC, sourceType);
            addPathConfig(configBuilder, PEACH_DIR_CFG, PEACH_DIR_DESC, sourceType);
            addPathConfig(configBuilder, VIRUS_DIR_CFG, VIRUS_DIR_DESC, sourceType);
            addPathConfig(configBuilder, CIDER_DIR_CFG, CIDER_DIR_DESC, sourceType);
            addPathConfig(configBuilder, TEAL_DIR_CFG, TEAL_DIR_DESC, sourceType);
            addPathConfig(configBuilder, V_CHORD_DIR_CFG, V_CHORD_DIR_DESC, sourceType);
            addPathConfig(configBuilder, SIGS_DIR_CFG, SIGS_DIR_DESC, sourceType);
            addPathConfig(configBuilder, ISOFOX_DIR_CFG, ISOFOX_DIR_DESC, sourceType);
            addPathConfig(configBuilder, TUMOR_FLAGSTAT, formSourceDescription("Tumor flagstat", sourceType), sourceType);
            addPathConfig(configBuilder, GERMLINE_FLAGSTAT, formSourceDescription("Germline flagstat", sourceType), sourceType);
            addPathConfig(configBuilder, TUMOR_BAM_METRICS, formSourceDescription("Tumor BAM metrics", sourceType), sourceType);
            addPathConfig(configBuilder, GERMLINE_BAM_METRICS, formSourceDescription("Germline BAM metrics", sourceType), sourceType);
            addPathConfig(configBuilder, SNP_GENOTYPE, formSourceDescription("SNP genotype", sourceType), sourceType);

            configBuilder.addPath(
                    formSourceConfig(SOMATIC_VCF, sourceType), false,
                    formSourceDescription("Non-Purple somatic VCF (eg Pave, Sage)", sourceType));

            configBuilder.addPath(
                    formSourceConfig(SOMATIC_UNFILTERED_VCF, sourceType), false,
                    formSourceDescription("VCF to search for filtered variants", sourceType));

            configBuilder.addConfigItem(
                    formSourceConfig(PIPELINE_FORMAT_CFG, sourceType), false,
                    formSourceDescription(PIPELINE_FORMAT_DESC, sourceType));
            configBuilder.addPath(
                    formSourceConfig(PIPELINE_FORMAT_FILE_CFG, sourceType), false,
                    formSourceDescription(PIPELINE_FORMAT_FILE_DESC, sourceType));
        }
    }

    private static String formSourceDescription(final String desc, final SourceType sourceType)
    {
        return format("%s: source %s", desc, sourceType.configStr());
    }

    private static String formSourceConfig(final String config, final SourceType sourceType)
    {
        return format("%s_%s", config, sourceType.configStr());
    }

    private static String getConfigValue(final ConfigBuilder configBuilder, final String config, SourceType sourceType)
    {
        return configBuilder.getValue(formSourceConfig(config, sourceType), "");
    }

    public static FileSources fromConfig(final SourceType sourceType, final ConfigBuilder configBuilder)
    {
        String sampleDir = checkAddDirSeparator(getConfigValue(configBuilder, SAMPLE_DIR, sourceType));

        PipelineToolDirectories defaultToolDirs = resolveDefaultToolDirs(configBuilder, sourceType);

        String linxDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.linxSomaticDir(), LINX_DIR_CFG, sourceType);

        String linxGermlineDir = getDirectory(
                configBuilder, sampleDir, defaultToolDirs.linxGermlineDir(), LINX_GERMLINE_DIR_CFG, sourceType);

        String cobaltDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.cobaltDir(), COBALT_DIR_CFG, sourceType);
        String purpleDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.purpleDir(), PURPLE_DIR_CFG, sourceType);
        String cuppaDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.cuppaDir(), CUPPA_DIR_CFG, sourceType);
        String lilacDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.lilacDir(), LILAC_DIR_CFG, sourceType);
        String chordDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.chordDir(), CHORD_DIR_CFG, sourceType);
        String peachDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.peachDir(), PEACH_DIR_CFG, sourceType);
        String virusDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.virusInterpreterDir(), VIRUS_DIR_CFG, sourceType);

        String somaticVcf = getConfigValue(configBuilder, SOMATIC_VCF, sourceType);
        String somaticUnfilteredVcf = getConfigValue(configBuilder, SOMATIC_UNFILTERED_VCF, sourceType);

        String tumorFlagstat = getDirectory(configBuilder, sampleDir, defaultToolDirs.tumorFlagstatDir(), TUMOR_FLAGSTAT, sourceType);
        String germlineFlagstat = getDirectory(configBuilder, sampleDir, defaultToolDirs.germlineFlagstatDir(), GERMLINE_FLAGSTAT, sourceType);
        String tumorBamMetrics = getDirectory(configBuilder, sampleDir, defaultToolDirs.tumorMetricsDir(), TUMOR_BAM_METRICS, sourceType);
        String germlineBamMetrics = getDirectory(configBuilder, sampleDir, defaultToolDirs.germlineMetricsDir(), GERMLINE_BAM_METRICS, sourceType);
        String snpGenotype = getDirectory(configBuilder, sampleDir, defaultToolDirs.snpGenotypeDir(), SNP_GENOTYPE, sourceType);
        String ciderDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.ciderDir(), CIDER_DIR_CFG, sourceType);
        String tealDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.tealDir(), TEAL_DIR_CFG, sourceType);
        String vChordDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.vChordDir(), V_CHORD_DIR_CFG, sourceType);
        String sigsDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.sigsDir(), SIGS_DIR_CFG, sourceType);
        String isofoxDir = getDirectory(configBuilder, sampleDir, defaultToolDirs.isofoxDir(), ISOFOX_DIR_CFG, sourceType);

        return new FileSources(sourceType, linxDir, cobaltDir, purpleDir, linxGermlineDir, cuppaDir, lilacDir, chordDir, peachDir, virusDir,
                somaticVcf, somaticUnfilteredVcf, tumorFlagstat, germlineFlagstat, tumorBamMetrics, germlineBamMetrics, snpGenotype,
                ciderDir, tealDir, vChordDir, sigsDir, isofoxDir);
    }

    private static PipelineToolDirectories resolveDefaultToolDirs(final ConfigBuilder configBuilder, final SourceType sourceType)
    {
        String pipelineFormatConfigStr = formSourceConfig(PIPELINE_FORMAT_CFG, sourceType);
        String pipelineFormatFileConfigStr = formSourceConfig(PIPELINE_FORMAT_FILE_CFG, sourceType);
        return PipelineToolDirectories.resolveToolDirectories(configBuilder, pipelineFormatConfigStr, pipelineFormatFileConfigStr);
    }

    private static String getDirectory(
            final ConfigBuilder configBuilder, final String sampleDir, final String toolDefaultDir,
            final String config, final SourceType sourceType)
    {
        // if a tool directory is specified in config, then it overrides the default pipeline directory
        // if the root sample directory is specified, then the tool directory is relative to that, otherwise is absolute

        String configStr = formSourceConfig(config, sourceType);

        if(!configBuilder.hasValue(configStr) && sampleDir.isEmpty())
            return "";

        String toolDir = configBuilder.getValue(configStr, toolDefaultDir);

        String directory = "";

        if(sampleDir.isEmpty())
            directory = toolDir;
        else if(toolDir.isEmpty())
            directory = sampleDir;
        else
            directory = format("%s%s", sampleDir, toolDir);

        return checkAddDirSeparator(directory);
    }
}

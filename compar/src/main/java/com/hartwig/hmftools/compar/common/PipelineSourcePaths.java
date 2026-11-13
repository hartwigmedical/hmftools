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
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PAVE_SOMATIC_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PAVE_SOMATIC_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PEACH_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAGE_DIR_DESC;
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

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.pipeline.PipelineToolDirectories;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class PipelineSourcePaths
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
    public final String SageSomatic;
    public final String PaveSomatic;
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
    private static final String SOMATIC_UNFILTERED_VCF = "somatic_unfiltered_vcf";
    private static final String TUMOR_FLAGSTAT = "tumor_flagstat_dir";
    private static final String GERMLINE_FLAGSTAT = "germline_flagstat_dir";
    private static final String TUMOR_BAM_METRICS = "tumor_bam_metrics_dir";
    private static final String GERMLINE_BAM_METRICS = "germline_bam_metrics_dir";
    private static final String SNP_GENOTYPE = "snp_genotype_dir";

    private final String mSampleDir;
    private final PipelineToolDirectories mDefaultToolDirs;
    private final Map<String,String> mSampleDirOverrides; // source sample ID to sample root directory

    private PipelineSourcePaths(
            final SourceType source, final String linx, final String cobalt, final String purple, final String linxGermline,
            final String cuppa, final String lilac, final String chord, final String peach, final String virus,
            final String sageSomaticDir, final String paveSomaticDir, final String somaticUnfilteredVcf,
            final String tumorFlagstat, final String germlineFlagstat, final String tumorBamMetrics,
            final String germlineBamMetrics, final String snpGenotype, final String cider, final String teal, final String vChord,
            final String sigs, final String isofox, final String sampleDir, final PipelineToolDirectories defaultToolDirs)
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
        SageSomatic = sageSomaticDir;
        PaveSomatic = paveSomaticDir;
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

        mSampleDir = sampleDir;
        mDefaultToolDirs = defaultToolDirs;
        mSampleDirOverrides = Maps.newHashMap();
    }

    public void addSampleDirOverride(final String sampleId, final String sampleDir)
    {
        mSampleDirOverrides.put(sampleId, sampleDir);
    }

    public static PipelineSourcePaths sampleInstance(final PipelineSourcePaths fileSources, final String sampleId, final String referenceId)
    {
        String sampleDir = checkAddDirSeparator(fileSources.mSampleDirOverrides.getOrDefault(sampleId, fileSources.mSampleDir));
        PipelineToolDirectories defaultToolDirs = fileSources.mDefaultToolDirs;

        return new PipelineSourcePaths(
                fileSources.Source,
                resolveSamplePath(sampleDir, defaultToolDirs.linxSomaticDir(), fileSources.Linx, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.cobaltDir(), fileSources.Cobalt, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.purpleDir(), fileSources.Purple, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.linxGermlineDir(), fileSources.LinxGermline, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.cuppaDir(), fileSources.Cuppa, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.lilacDir(), fileSources.Lilac, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.chordDir(), fileSources.Chord, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.peachDir(), fileSources.Peach, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.virusInterpreterDir(), fileSources.Virus, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.sageSomaticDir(), fileSources.SageSomatic, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.paveSomaticDir(), fileSources.PaveSomatic, sampleId, referenceId),
                convertWildcardSamplePath(fileSources.SomaticUnfilteredVcf, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.tumorFlagstatDir(), fileSources.TumorFlagstat, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.germlineFlagstatDir(), fileSources.GermlineFlagstat, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.tumorMetricsDir(), fileSources.TumorBamMetrics, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.germlineMetricsDir(), fileSources.GermlineBamMetrics, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.snpGenotypeDir(), fileSources.SnpGenotype, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.ciderDir(), fileSources.Cider, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.tealDir(), fileSources.Teal, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.vChordDir(), fileSources.VChord, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.sigsDir(), fileSources.Sigs, sampleId, referenceId),
                resolveSamplePath(sampleDir, defaultToolDirs.isofoxDir(), fileSources.Isofox, sampleId, referenceId),
                null, null);
    }

    private static String resolveSamplePath(
            final String sampleDir, final String toolDefaultDir, final String configuredToolDir,
            final String sampleId, final String referenceId)
    {
        return convertWildcardSamplePath(getDirectory(sampleDir, toolDefaultDir, configuredToolDir), sampleId, referenceId);
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
            addPathConfig(configBuilder, SAGE_DIR_CFG, SAGE_DIR_DESC, sourceType);
            addPathConfig(configBuilder, PAVE_SOMATIC_DIR_CFG, PAVE_SOMATIC_DIR_DESC, sourceType);
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

    public static PipelineSourcePaths fromConfig(final SourceType sourceType, final ConfigBuilder configBuilder)
    {
        String sampleDir = getConfigValue(configBuilder, SAMPLE_DIR, sourceType);

        PipelineToolDirectories defaultToolDirs = resolveDefaultToolDirs(configBuilder, sourceType);

        return new PipelineSourcePaths(
                sourceType,
                getConfiguredToolDir(configBuilder, LINX_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, COBALT_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, PURPLE_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, LINX_GERMLINE_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, CUPPA_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, LILAC_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, CHORD_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, PEACH_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, VIRUS_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, SAGE_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, PAVE_SOMATIC_DIR_CFG, sourceType),
                getConfigValue(configBuilder, SOMATIC_UNFILTERED_VCF, sourceType),
                getConfiguredToolDir(configBuilder, TUMOR_FLAGSTAT, sourceType),
                getConfiguredToolDir(configBuilder, GERMLINE_FLAGSTAT, sourceType),
                getConfiguredToolDir(configBuilder, TUMOR_BAM_METRICS, sourceType),
                getConfiguredToolDir(configBuilder, GERMLINE_BAM_METRICS, sourceType),
                getConfiguredToolDir(configBuilder, SNP_GENOTYPE, sourceType),
                getConfiguredToolDir(configBuilder, CIDER_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, TEAL_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, V_CHORD_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, SIGS_DIR_CFG, sourceType),
                getConfiguredToolDir(configBuilder, ISOFOX_DIR_CFG, sourceType),
                sampleDir, defaultToolDirs);
    }

    private static String getConfiguredToolDir(final ConfigBuilder configBuilder, final String config, final SourceType sourceType)
    {
        String configStr = formSourceConfig(config, sourceType);
        return configBuilder.hasValue(configStr) ? configBuilder.getValue(configStr) : null;
    }

    private static PipelineToolDirectories resolveDefaultToolDirs(final ConfigBuilder configBuilder, final SourceType sourceType)
    {
        String pipelineFormatConfigStr = formSourceConfig(PIPELINE_FORMAT_CFG, sourceType);
        String pipelineFormatFileConfigStr = formSourceConfig(PIPELINE_FORMAT_FILE_CFG, sourceType);
        return PipelineToolDirectories.resolveToolDirectories(configBuilder, pipelineFormatConfigStr, pipelineFormatFileConfigStr);
    }

    private static String getDirectory(final String sampleDir, final String toolDefaultDir, final String configuredToolDir)
    {
        // if a tool directory is specified in config, then it overrides the default pipeline directory
        // if the root sample directory is specified, then the tool directory is relative to that, otherwise is absolute

        if(configuredToolDir == null && sampleDir.isEmpty())
            return "";

        String toolDir = configuredToolDir != null ? configuredToolDir : toolDefaultDir;

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

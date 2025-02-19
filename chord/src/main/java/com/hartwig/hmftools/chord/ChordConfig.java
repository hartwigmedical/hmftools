package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SOMATIC_VCF_SUFFIX;
import static com.hartwig.hmftools.common.purple.PurpleCommon.PURPLE_SV_VCF_SUFFIX;
import static com.hartwig.hmftools.common.utils.TaskExecutor.THREADS;
import static com.hartwig.hmftools.common.utils.TaskExecutor.THREADS_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.TaskExecutor;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;


public class ChordConfig
{
    public final List<String> SampleIds;

    public final String PurpleDir;
    public final String SnvIndelVcfFile;
    public final String SvVcfFile;

    public final String RefGenomeFile;

    public final String OutputDir;
    public final String OutputId;

    public final int Threads;

    public final boolean IncludeNonPass;
    public final boolean WriteDetailedFiles;

    private static final String SNV_INDEL_VCF_FILE = "snv_indel_vcf_file";
    private static final String SNV_INDEL_VCF_FILE_DESC = "Path to the VCF containing SNVs and indels";

    private static final String SV_VCF_FILE = "sv_vcf_file";
    private static final String SV_VCF_FILE_DESC = "Path to the VCF containing structural variants";

    private static final String INCLUDE_NON_PASS = "include_non_pass";
    private static final String INCLUDE_NON_PASS_DESC = "Include non pass variants when counting mutation types";

    private static final String WRITE_DETAILED_FILES = "write_detailed_files";
    private static final String WRITE_DETAILED_FILES_DESC = "Write mutation context attributes per variant";

    private static final String SAMPLE_IDS_DELIM = ",";

    public ChordConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = getSampleIds(configBuilder);

        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        SnvIndelVcfFile = configBuilder.getValue(SNV_INDEL_VCF_FILE);
        SvVcfFile = configBuilder.getValue(SV_VCF_FILE);
        checkRequiredInputPaths(PurpleDir, SnvIndelVcfFile, SvVcfFile);

        RefGenomeFile = configBuilder.getValue(REF_GENOME);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        Threads = TaskExecutor.parseThreads(configBuilder);

        IncludeNonPass = configBuilder.hasFlag(INCLUDE_NON_PASS);
        WriteDetailedFiles = configBuilder.hasFlag(WRITE_DETAILED_FILES);
    }

    private List<String> getSampleIds(ConfigBuilder configBuilder)
    {
        if(configBuilder.hasValue(SAMPLE))
        {
            String sampleIdsString = configBuilder.getValue(SAMPLE);
            return List.of(sampleIdsString.split(SAMPLE_IDS_DELIM));
        }

        if(configBuilder.hasValue(SAMPLE_ID_FILE))
        {
            return ConfigUtils.loadSampleIdsFile(configBuilder);
        }

        CHORD_LOGGER.error("Either -{} or -{} must be provided to config", SAMPLE, SAMPLE_ID_FILE);
        System.exit(1);
        return null;
    }

    public boolean isMultiSample(){ return SampleIds.size() > 1; }

    public boolean isSingleSample(){ return !isMultiSample(); }

    public String snvIndelVcfFile(final String sampleId)
    {
        return (SnvIndelVcfFile != null) ?
                ConfigUtils.convertWildcardSamplePath(SnvIndelVcfFile, sampleId) :
                ConfigUtils.convertWildcardSamplePath(PurpleDir + "/*" + PURPLE_SOMATIC_VCF_SUFFIX, sampleId);
    }

    public String svVcfFile(final String sampleId)
    {
        return (SvVcfFile != null) ?
                ConfigUtils.convertWildcardSamplePath(SvVcfFile, sampleId) :
                ConfigUtils.convertWildcardSamplePath(PurpleDir + "/*" + PURPLE_SV_VCF_SUFFIX, sampleId);
    }

    public static void checkRequiredInputPaths(String purpleDir, String snvIndelVcfFile, String svVcfFile)
    {
        if(purpleDir != null)
            return;

        if(snvIndelVcfFile != null && svVcfFile != null)
            return;

        CHORD_LOGGER.error("Please provide either 1) {}, or 2) {} and {}", PURPLE_DIR_CFG, SNV_INDEL_VCF_FILE, SV_VCF_FILE);
        System.exit(1);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);

        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(SNV_INDEL_VCF_FILE, false, SNV_INDEL_VCF_FILE_DESC);
        configBuilder.addPath(SV_VCF_FILE, false, SV_VCF_FILE_DESC);

        configBuilder.addConfigItem(REF_GENOME, false, REF_GENOME_CFG_DESC);

        FileWriterUtils.addOutputOptions(configBuilder);

        configBuilder.addConfigItem(THREADS, false, THREADS_DESC, "1");

        configBuilder.addFlag(INCLUDE_NON_PASS, INCLUDE_NON_PASS_DESC);
        configBuilder.addFlag(WRITE_DETAILED_FILES, WRITE_DETAILED_FILES_DESC);

        ConfigUtils.addLoggingOptions(configBuilder);
    }

    @VisibleForTesting
    public ChordConfig(
            List<String> sampleIds, String purpleDir, String snvIndelVcfFile, String svFileFile, String refGenomeFile,
            String outputDir, String outputId, int threads, boolean includeNonPass, boolean writeDetailedFiles
    )
    {
        SampleIds = sampleIds;
        PurpleDir = purpleDir;
        SnvIndelVcfFile = snvIndelVcfFile;
        SvVcfFile = svFileFile;
        RefGenomeFile = refGenomeFile;
        OutputDir = outputDir;
        OutputId = outputId;
        Threads = threads;
        IncludeNonPass = includeNonPass;
        WriteDetailedFiles = writeDetailedFiles;
    }

    @VisibleForTesting
    public static class Builder
    {
        private List<String> SampleIds;
        private String PurpleDir;
        private String SnvIndelVcfFile;
        private String SvVcfFile;
        private String RefGenomeFile;
        private String OutputDir;
        private String OutputId = "";
        private int Threads = 1;
        private boolean IncludeNonPass = false;
        private boolean WriteDetailedFiles = false;

        public Builder sampleIds(String sampleId)
        {
            SampleIds = List.of(sampleId);
            return this;
        }

        public Builder sampleIds(List<String> sampleIds)
        {
            SampleIds = sampleIds;
            return this;
        }

        public Builder purpleDir(String purpleDir)
        {
            PurpleDir = purpleDir;
            return this;
        }

        public Builder snvIndelVcfFile(String snvIndelVcfFile)
        {
            SnvIndelVcfFile = snvIndelVcfFile;
            return this;
        }

        public Builder svVcfFile(String svVcfFile)
        {
            SvVcfFile = svVcfFile;
            return this;
        }

        public Builder refGenomeFile(String refGenomeFile)
        {
            RefGenomeFile = refGenomeFile;
            return this;
        }

        public Builder outputDir(String outputDir)
        {
            OutputDir = outputDir;
            return this;
        }

        public Builder outputId(String outputId)
        {
            OutputId = outputId;
            return this;
        }

        public Builder threads(int threads)
        {
            Threads = threads;
            return this;
        }

        public Builder includeNonPass(boolean includeNonPass)
        {
            IncludeNonPass = includeNonPass;
            return this;
        }

        public Builder writeDetailedFiles(boolean writeDetailedFiles)
        {
            WriteDetailedFiles = writeDetailedFiles;
            return this;
        }

        public ChordConfig build()
        {
            return new ChordConfig(
                    SampleIds, PurpleDir, SnvIndelVcfFile, SvVcfFile, RefGenomeFile,
                    OutputDir, OutputId, Threads, IncludeNonPass, WriteDetailedFiles
            );
        }
    }
}

package com.hartwig.hmftools.chord;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;


public class ChordConfig
{
    public final List<String> SampleIds;

    public final String PurpleDir;

    public final RefGenomeVersion RefGenVersion;

    public final String OutputDir;
    public final String OutputId;

    public final boolean IncludeNonPass;

    public final String ChordToolDir;

    private static final String INCLUDE_NON_PASS = "include_non_pass";
    private static final String INCLUDE_NON_PASS_DESC = "Include non pass variants when counting mutation types";

    private static final String CHORD_TOOL_DIR = "chord_tool_dir";
    private static final String CHORD_TOOL_DIR_DESC = "Dir containing the CHORD and mutSigExtractor R packages";

    private static final String SAMPLE_IDS_DELIM = ",";

    public ChordConfig(final ConfigBuilder configBuilder)
    {
        SampleIds = parseSampleIds(configBuilder.getValue(SAMPLE));

        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        IncludeNonPass = configBuilder.hasFlag(INCLUDE_NON_PASS);

        ChordToolDir = configBuilder.getValue(CHORD_TOOL_DIR);
    }

    private List<String> parseSampleIds(String sampleIdsString)
    {
        List<String> sampleIds = List.of(sampleIdsString.split(SAMPLE_IDS_DELIM));

        if(sampleIds.size() > 1)
        {
            throw new UnsupportedOperationException("Running CHORD in multi-sample mode is not yet implemented");
        }

        return sampleIds;
    }

    public String purpleSomaticVcfFile(String sampleId) { return PurpleCommon.purpleSomaticVcfFile(PurpleDir, sampleId); }

    public String purpleSvVcfFile(String sampleId) { return PurpleCommon.purpleSomaticSvFile(PurpleDir, sampleId); }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);

        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);

        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());

        configBuilder.addFlag(INCLUDE_NON_PASS, INCLUDE_NON_PASS_DESC);

        configBuilder.addPath(CHORD_TOOL_DIR, false, CHORD_TOOL_DIR_DESC);

        FileWriterUtils.addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }

    public ChordConfig(
            List<String> sampleIds, String purpleDir, RefGenomeVersion refGenomeVersion, String outputDir, String outputId,
            boolean includeNonPass, String chordToolDir)
    {
        SampleIds = sampleIds;
        PurpleDir = purpleDir;
        RefGenVersion = refGenomeVersion;
        OutputDir = outputDir;
        OutputId = outputId;
        IncludeNonPass = includeNonPass;
        ChordToolDir = chordToolDir;
    }

    public static class Builder
    {
        private List<String> SampleIds;
        private String PurpleDir;
        private RefGenomeVersion RefGenVersion = V37;
        private String OutputDir;
        private String OutputId = "";
        private boolean IncludeNonPass = false;
        private String ChordToolDir = "";

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

        public Builder refGenomeVersion(RefGenomeVersion refGenomeVersion)
        {
            RefGenVersion = refGenomeVersion;
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

        public Builder includeNonPass(boolean includeNonPass)
        {
            IncludeNonPass = includeNonPass;
            return this;
        }

        public Builder chordToolDir(String chordToolDir)
        {
            ChordToolDir = chordToolDir;
            return this;
        }

        public ChordConfig build(){ return new ChordConfig(SampleIds, PurpleDir, RefGenVersion, OutputDir, OutputId, IncludeNonPass, ChordToolDir); }
    }
}

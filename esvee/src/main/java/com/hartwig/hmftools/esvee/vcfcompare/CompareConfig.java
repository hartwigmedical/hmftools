package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigItem.enumValueSelectionAsStr;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.esvee.common.FileCommon.ESVEE_FILE_ID;
import static com.hartwig.hmftools.esvee.common.FileCommon.FILE_NAME_DELIM;
import static com.hartwig.hmftools.esvee.common.FileCommon.formOutputFile;
import static com.hartwig.hmftools.esvee.vcfcompare.CompareTask.MATCH_BREAKENDS;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.common.utils.file.FileWriterUtils;
import com.hartwig.hmftools.esvee.assembly.output.WriteType;

import org.jetbrains.annotations.Nullable;

public class CompareConfig
{
    public final List<CompareTask> CompareTasks;

    public final String SampleId;
    public final String ReferenceId;

    public final String OldVcf;
    public final String NewVcf;

    public final String OldUnfilteredVcf;
    public final String NewUnfilteredVcf;

    public final String OutputDir;
    public final String OutputId;

    public final boolean IncludeNonPass;
    public final boolean WriteMatches;
    public RefGenomeVersion RefGenVersion;

    static final String OLD_VCF = "old_vcf";
    static final String NEW_VCF = "new_vcf";

    static final String OLD_UNFILTERED_VCF = "old_unfiltered_vcf";
    static final String NEW_UNFILTERED_VCF = "new_unfiltered_vcf";

    private static final String INCLUDE_NON_PASS = "include_non_pass";
    private static final String WRITE_MATCHES = "write_matches";

    private static final String COMPARE_TASKS = "tasks";

    public CompareConfig(final ConfigBuilder configBuilder)
    {
        CompareTasks = CompareTask.fromConfig(configBuilder.getValue(COMPARE_TASKS));

        SampleId = configBuilder.getValue(SAMPLE);
        ReferenceId = configBuilder.getValue(REFERENCE, "");

        OldVcf = configBuilder.getValue(OLD_VCF);
        NewVcf = configBuilder.getValue(NEW_VCF);

        OldUnfilteredVcf = configBuilder.getValue(OLD_UNFILTERED_VCF);
        NewUnfilteredVcf = configBuilder.getValue(NEW_UNFILTERED_VCF);

        OutputDir = FileWriterUtils.parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        IncludeNonPass = configBuilder.hasFlag(INCLUDE_NON_PASS);
        WriteMatches = configBuilder.hasFlag(WRITE_MATCHES);
        RefGenVersion = RefGenomeVersion.from(configBuilder);
    }

    public String formFilename(final String fileType)
    {
        String appStage = ESVEE_FILE_ID + FILE_NAME_DELIM + "compare" + FILE_NAME_DELIM + fileType;
        return formOutputFile(OutputDir, SampleId, appStage, "tsv", OutputId);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(
                COMPARE_TASKS, false, enumValueSelectionAsStr(WriteType.values(), "Comparison types"),
                MATCH_BREAKENDS.toString());

        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addConfigItem(REFERENCE, false, REFERENCE_DESC);

        configBuilder.addPath(OLD_VCF, true, "Path to the old VCF file");
        configBuilder.addPath(NEW_VCF, true, "Path to the new VCF file");

        configBuilder.addPath(OLD_UNFILTERED_VCF, false, "Path to the old unfiltered VCF file");
        configBuilder.addPath(NEW_UNFILTERED_VCF, false, "Path to the new unfiltered VCF file");

        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());
        configBuilder.addFlag(INCLUDE_NON_PASS, "Show variants not PASSing in both old nor new VCF files");
        configBuilder.addFlag(WRITE_MATCHES, "Write exact matches to output");

        FileWriterUtils.addOutputOptions(configBuilder);
        ConfigUtils.addLoggingOptions(configBuilder);
    }
}

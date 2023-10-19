package com.hartwig.hmftools.cup.prep;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.ISOFOX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.LINX_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.VIRUS_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.SAMPLE_ID_FILE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.cup.CuppaConfig.CATEGORIES;
import static com.hartwig.hmftools.cup.CuppaConfig.configCategories;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class PrepConfig
{
    public final RefGenomeVersion RefGenVersion;
    public final List<CategoryType> Categories;

    public final List<String> SampleIds;

    public final String OutputDir;
    public final String OutputId; // for multi-sample mode
    public final boolean WriteByCategory;

    // pipeline directories, accepting wildcards
    public final String LinxDir;
    public final String PurpleDir;
    public final String VirusDir;
    public final String IsofoxDir;

    private static final String WRITE_FILE_BY_CATEGORY = "write_by_category";

    public PrepConfig(final ConfigBuilder configBuilder)
    {
        Categories = configCategories(configBuilder);

        SampleIds = Lists.newArrayList();

        if(configBuilder.hasValue(SAMPLE))
        {
            SampleIds.add(configBuilder.getValue(SAMPLE));
        }
        else
        {
            SampleIds.addAll(loadSampleIdsFile(configBuilder));
        }

        RefGenVersion = RefGenomeVersion.from(configBuilder);

        LinxDir = configBuilder.getValue(LINX_DIR_CFG, "");
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG, "");
        VirusDir = configBuilder.getValue(VIRUS_DIR_CFG, "");
        IsofoxDir = configBuilder.getValue(ISOFOX_DIR_CFG, "");

        OutputDir = parseOutputDir(configBuilder);
        OutputId = configBuilder.getValue(OUTPUT_ID);

        WriteByCategory = SampleIds.size() > 1 && configBuilder.hasFlag(WRITE_FILE_BY_CATEGORY);
    }

    public boolean isMultiSample() { return SampleIds.size() > 1; }
    public boolean isSingleSample() { return SampleIds.size() == 1; }

    public static void addPipelineDirectories(final ConfigBuilder configBuilder)
    {
        configBuilder.addPath(LINX_DIR_CFG, false, LINX_DIR_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(VIRUS_DIR_CFG, false, PURPLE_DIR_DESC);
        configBuilder.addPath(ISOFOX_DIR_CFG, false, ISOFOX_DIR_DESC);
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(CATEGORIES, false, "Categories to build ref data for");

        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        configBuilder.addPath(SAMPLE_ID_FILE, false, SAMPLE_ID_FILE_DESC);
        configBuilder.addConfigItem(REF_GENOME_VERSION, false, REF_GENOME_VERSION_CFG_DESC, V37.toString());
        configBuilder.addFlag(WRITE_FILE_BY_CATEGORY, "Cohort mode - write files by category");

        addPipelineDirectories(configBuilder);

        /*
        configBuilder.addConfigItem(REF_SNV_SAMPLE_POS_FREQ_FILE, false, "Ref SNV position frequency matrix data file");
        configBuilder.addConfigItem(REF_SNV_COUNTS_FILE, false, "Ref SNV trinucleotide matrix data file");
        configBuilder.addPath(REF_GENE_EXP_DATA_FILE, false, "Ref sample RNA gene expression cohort data file");
        configBuilder.addPath(REF_ALT_SJ_DATA_FILE, false, "Ref sample RNA alternate splice junction cohort data file");
        */

        addLoggingOptions(configBuilder);
        addOutputOptions(configBuilder);
    }
}

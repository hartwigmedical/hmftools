package com.hartwig.hmftools.virusinterpreter;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_CFG;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.PURPLE_DIR_DESC;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.addOutputDir;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.parseOutputDir;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class VirusInterpreterConfig
{
    public final String SampleId;
    public final String VirusBreakendTsv;
    public final String TaxonomyDbTsv;
    public final String VirusReportedDbTsv;
    public final String VirusBlacklistedDbTsv;
    public final String PurpleDir;
    public final String TumorSampleWGSMetricsFile;
    public final String OutputDir;

    private static final String VIRUS_BREAKEND_TSV = "virus_breakend_tsv";
    private static final String TAXONOMY_DB_TSV = "taxonomy_db_tsv";
    private static final String VIRUS_REPORTING_DB_TSV = "virus_reporting_db_tsv";
    private static final String VIRUS_BLACKLISTING_DB_TSV = "virus_blacklisting_db_tsv";
    private static final String TUMOR_SAMPLE_WGS_METRICS_FILE = "tumor_sample_wgs_metrics_file";

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addConfigItem(SAMPLE, true, SAMPLE_DESC);
        configBuilder.addPath(PURPLE_DIR_CFG, true, PURPLE_DIR_DESC);
        configBuilder.addPath(VIRUS_BREAKEND_TSV, true, "Path towards the virus breakend TSV");
        configBuilder.addPath(TAXONOMY_DB_TSV, true, "Path towards a TSV containing a mapping from taxid to taxonomy name");
        configBuilder.addPath(VIRUS_REPORTING_DB_TSV, true, "Path towards a TSV containing reported of viruses");
        configBuilder.addPath(VIRUS_BLACKLISTING_DB_TSV, true, "Path towards a TSV containing viruses to blacklist");
        configBuilder.addPath(TUMOR_SAMPLE_WGS_METRICS_FILE, true, "Path towards the tumor sample WGS metrics file");

        addOutputDir(configBuilder);
        addLoggingOptions(configBuilder);
    }

    public VirusInterpreterConfig(final ConfigBuilder configBuilder)
    {
        SampleId = configBuilder.getValue(SAMPLE);
        PurpleDir = configBuilder.getValue(PURPLE_DIR_CFG);
        VirusBreakendTsv = configBuilder.getValue(VIRUS_BREAKEND_TSV);
        TaxonomyDbTsv = configBuilder.getValue(TAXONOMY_DB_TSV);
        VirusReportedDbTsv = configBuilder.getValue(VIRUS_REPORTING_DB_TSV);
        VirusBlacklistedDbTsv = configBuilder.getValue(VIRUS_BLACKLISTING_DB_TSV);
        TumorSampleWGSMetricsFile = configBuilder.getValue(TUMOR_SAMPLE_WGS_METRICS_FILE);
        OutputDir = parseOutputDir(configBuilder);
    }
}

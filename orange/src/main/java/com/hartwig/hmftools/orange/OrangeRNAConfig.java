package com.hartwig.hmftools.orange;

import com.hartwig.hmftools.orange.util.Config;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface OrangeRNAConfig {

    Logger LOGGER = LogManager.getLogger(OrangeRNAConfig.class);

    String RNA_SAMPLE_ID = "rna_sample_id";

    String ISOFOX_GENE_DISTRIBUTION_CSV = "isofox_gene_distribution_csv";
    String ISOFOX_ALT_SJ_COHORT_CSV = "isofox_alt_sj_cohort_csv";

    String ISOFOX_SUMMARY_CSV = "isofox_summary_csv";
    String ISOFOX_GENE_DATA_CSV = "isofox_gene_data_csv";
    String ISOFOX_FUSION_CSV = "isofox_fusion_csv";
    String ISOFOX_ALT_SPLICE_JUNCTION_CSV = "isofox_alt_splice_junction_csv";

    @NotNull
    static Options createOptions() {
        Options options = new Options();

        options.addOption(RNA_SAMPLE_ID, true, "(Optional) The RNA sample of the tumor sample for which ORANGE will run.");

        options.addOption(ISOFOX_GENE_DISTRIBUTION_CSV, true, "Path to isofox gene distribution CSV.");
        options.addOption(ISOFOX_ALT_SJ_COHORT_CSV, true, "Path to isofox alt SJ cohort CSV.");

        options.addOption(ISOFOX_SUMMARY_CSV, true, "Path towards the ISOFOX summary data.");
        options.addOption(ISOFOX_GENE_DATA_CSV, true, "Path towards the ISOFOX gene data.");
        options.addOption(ISOFOX_FUSION_CSV, true, "Path towards the ISOFOX fusion data.");
        options.addOption(ISOFOX_ALT_SPLICE_JUNCTION_CSV, true, "Path towards the ISOFOX alt splice junction data.");

        return options;
    }

    @NotNull
    String rnaSampleId();

    @NotNull
    String isofoxGeneDistributionCsv();

    @NotNull
    String isofoxAltSjCohortCsv();

    @NotNull
    String isofoxSummaryCsv();

    @NotNull
    String isofoxGeneDataCsv();

    @NotNull
    String isofoxFusionCsv();

    @NotNull
    String isofoxAltSpliceJunctionCsv();

    @NotNull
    static OrangeRNAConfig createConfig(@NotNull CommandLine cmd) throws ParseException {
        String rnaSampleId = Config.optionalValue(cmd, RNA_SAMPLE_ID);
        if (rnaSampleId != null) {
            LOGGER.debug("RNA sample configured to {}", rnaSampleId);
        }

        if (rnaSampleId == null) {
            return null;
        }

        return ImmutableOrangeRNAConfig.builder()
                .rnaSampleId(rnaSampleId)
                .isofoxGeneDistributionCsv(Config.nonOptionalFile(cmd, ISOFOX_GENE_DISTRIBUTION_CSV))
                .isofoxAltSjCohortCsv(Config.nonOptionalFile(cmd, ISOFOX_ALT_SJ_COHORT_CSV))
                .isofoxSummaryCsv(Config.nonOptionalFile(cmd, ISOFOX_SUMMARY_CSV))
                .isofoxGeneDataCsv(Config.nonOptionalFile(cmd, ISOFOX_GENE_DATA_CSV))
                .isofoxFusionCsv(Config.nonOptionalFile(cmd, ISOFOX_FUSION_CSV))
                .isofoxAltSpliceJunctionCsv(Config.nonOptionalFile(cmd, ISOFOX_ALT_SPLICE_JUNCTION_CSV))
                .build();
    }
}

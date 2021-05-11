package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.patientreporter.germline.GermlineCondition;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingEntry;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.germline.ImmutableGermlineReportingEntry;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlackListModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlacklistFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryfile;

import org.junit.Test;

public class GenomicAnalyzerTest {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String PURPLE_PURITY_TSV = BASE_DIRECTORY + "/purple/sample.purple.purity.tsv";
    private static final String PURPLE_QC_FILE = BASE_DIRECTORY + "/purple/sample.purple.qc";
    private static final String PURPLE_DRIVER_CATALOG_SOMATIC_TSV = BASE_DIRECTORY + "/purple/sample.driver.catalog.somatic.tsv";
    private static final String PURPLE_DRIVER_CATALOG_GERMLINE_TSV = BASE_DIRECTORY + "/purple/sample.driver.catalog.germline.tsv";
    private static final String PURPLE_SOMATIC_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.somatic.vcf";
    private static final String PURPLE_GERMLINE_VARIANT_VCF = BASE_DIRECTORY + "/purple/sample.purple.germline.vcf";
    private static final String LINX_FUSIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.fusion.tsv";
    private static final String LINX_BREAKEND_TSV = BASE_DIRECTORY + "/linx/sample.linx.breakend.tsv";
    private static final String LINX_DRIVERS_TSV = BASE_DIRECTORY + "/linx/sample.linx.driver.catalog.tsv";
    private static final String CHORD_PREDICTION_TXT = BASE_DIRECTORY + "/chord/sample_chord_prediction.txt";
    private static final String PROTECT_EVIDENCE_TSV = BASE_DIRECTORY + "/protect/sample.protect.tsv";
    private static final String VIRUS_BREAKEND_TSV = BASE_DIRECTORY + "/virusbreakend/sample.virusbreakend.vcf.summary.tsv";
    private static final String PEACH_GENOTYPE_TSV = BASE_DIRECTORY + "/peach/sample.peach.genotype.tsv";
    private static final String VIRUS_DB_TSV = Resources.getResource("virusbreakend/virusdb.tsv").getPath();
    private static final String VIRUS_SUMMARY_TSV = Resources.getResource("virusbreakend/virus_summary.tsv").getPath();
    private static final String VIRUS_BLACKLIST_TSV = Resources.getResource("virusbreakend/virus_blacklist.tsv").getPath();

    @Test
    public void canRunOnTestRun() throws IOException {
        GenomicAnalyzer analyzer = new GenomicAnalyzer();

        GermlineReportingEntry germlineReportingTrue = ImmutableGermlineReportingEntry.builder()
                .gene("GENE1")
                .notifyClinicalGeneticist(GermlineCondition.ALWAYS)
                .conditionFilter(null)
                .build();

        GermlineReportingEntry germlineReportingFalse = ImmutableGermlineReportingEntry.builder()
                .gene("GENE2")
                .notifyClinicalGeneticist(GermlineCondition.NEVER)
                .conditionFilter(null)
                .build();

        GermlineReportingModel victim = new GermlineReportingModel(Lists.newArrayList(germlineReportingTrue, germlineReportingFalse));

        VirusDbModel virusDbModel = VirusDbFile.buildFromTsv(VIRUS_DB_TSV);
        VirusSummaryModel virusSummaryModel = VirusSummaryfile.buildFromTsv(VIRUS_SUMMARY_TSV);

        VirusBlackListModel virusBlackListModel = VirusBlacklistFile.buildFromTsv(VIRUS_BLACKLIST_TSV);

        assertNotNull(analyzer.run("sample",
                PURPLE_PURITY_TSV,
                PURPLE_QC_FILE,
                PURPLE_DRIVER_CATALOG_SOMATIC_TSV,
                PURPLE_DRIVER_CATALOG_GERMLINE_TSV,
                PURPLE_SOMATIC_VARIANT_VCF,
                PURPLE_GERMLINE_VARIANT_VCF,
                LINX_FUSIONS_TSV,
                LINX_BREAKEND_TSV,
                LINX_DRIVERS_TSV,
                CHORD_PREDICTION_TXT,
                PROTECT_EVIDENCE_TSV,
                VIRUS_BREAKEND_TSV,
                PEACH_GENOTYPE_TSV,
                victim,
                LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                virusDbModel,
                virusSummaryModel, virusBlackListModel));
    }
}

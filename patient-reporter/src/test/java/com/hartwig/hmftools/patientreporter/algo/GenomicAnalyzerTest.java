package com.hartwig.hmftools.patientreporter.algo;

import static org.junit.Assert.assertNotNull;

import java.io.IOException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.patientreporter.PatientReporterConfig;
import com.hartwig.hmftools.patientreporter.PatientReporterTestConfig;
import com.hartwig.hmftools.patientreporter.germline.GermlineCondition;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingEntry;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.germline.ImmutableGermlineReportingEntry;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlackListModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlacklistFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryModel;

import org.junit.Test;

public class GenomicAnalyzerTest {

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

        GermlineReportingModel germlineReportingModel =
                new GermlineReportingModel(Lists.newArrayList(germlineReportingTrue, germlineReportingFalse));

        PatientReporterConfig config = PatientReporterTestConfig.create();
        VirusDbModel virusDbModel = VirusDbFile.buildFromTsv(config.virusDbTsv());
        VirusSummaryModel virusSummaryModel = VirusSummaryFile.buildFromTsv(config.virusSummaryTsv());
        VirusBlackListModel virusBlackListModel = VirusBlacklistFile.buildFromTsv(config.virusBlacklistTsv());

        assertNotNull(analyzer.run("sample",
                config,
                germlineReportingModel,
                LimsGermlineReportingLevel.REPORT_WITH_NOTIFICATION,
                virusDbModel,
                virusSummaryModel,
                virusBlackListModel));
    }
}

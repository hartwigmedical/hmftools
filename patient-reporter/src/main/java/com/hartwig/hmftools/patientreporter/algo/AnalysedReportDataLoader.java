package com.hartwig.hmftools.patientreporter.algo;

import java.io.IOException;

import com.hartwig.hmftools.patientreporter.ReportData;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingFile;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.patientreporter.summary.SummaryFile;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusBlacklistFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusDbFile;
import com.hartwig.hmftools.patientreporter.virusbreakend.VirusSummaryFile;

import org.jetbrains.annotations.NotNull;

public final class AnalysedReportDataLoader {

    private AnalysedReportDataLoader() {
    }

    @NotNull
    public static AnalysedReportData buildFromFiles(@NotNull ReportData reportData, @NotNull String germlineReportingTsv,
            @NotNull String sampleSummaryTsv, @NotNull String virusTsv, @NotNull String virusSummaryTsv, @NotNull String virusBlacklistTsv)
            throws IOException {
        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromTsv(germlineReportingTsv);
        SummaryModel summaryModel = SummaryFile.buildFromTsv(sampleSummaryTsv);

        return ImmutableAnalysedReportData.builder()
                .from(reportData)
                .germlineReportingModel(germlineReportingModel)
                .summaryModel(summaryModel)
                .virusDbModel(VirusDbFile.buildFromTsv(virusTsv))
                .virusSummaryModel(VirusSummaryFile.buildFromTsv(virusSummaryTsv))
                .virusBlackListModel(VirusBlacklistFile.buildFromTsv(virusBlacklistTsv))
                .build();
    }
}

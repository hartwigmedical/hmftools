package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneFile;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelAssembly;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelFactory;
import com.hartwig.hmftools.patientreporter.summary.SummaryFile;
import com.hartwig.hmftools.patientreporter.summary.SummaryModel;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingFile;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineReportingModel;

import org.jetbrains.annotations.NotNull;

final class AnalysedReportDataLoader {

    private AnalysedReportDataLoader() {
    }

    @NotNull
    static AnalysedReportData buildFromFiles(@NotNull ReportData reportData, @NotNull String knowledgebaseDir,
            @NotNull String germlineGenesCsv, @NotNull String sampleSummaryTsv, @NotNull String driverGenePanelFile) throws IOException {
        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(knowledgebaseDir);

        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(germlineGenesCsv);
        SummaryModel summaryModel = SummaryFile.buildFromTsv(sampleSummaryTsv);
        List<DriverGene> driverGenes = DriverGeneFile.read(driverGenePanelFile);

        return ImmutableAnalysedReportData.builder()
                .from(reportData)
                .driverGene(driverGenes)
                .actionabilityAnalyzer(actionabilityAnalyzer)
                .germlineReportingModel(germlineReportingModel)
                .summaryModel(summaryModel)
                .build();
    }
}

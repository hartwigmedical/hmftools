package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;

import java.io.IOException;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.patientreporter.SequencedPatientReport;
import com.hartwig.hmftools.patientreporter.report.data.EvidenceReportData;
import com.hartwig.hmftools.patientreporter.report.data.ImmutableEvidenceReportData;
import com.hartwig.hmftools.patientreporter.report.pages.AlterationDebugPage;
import com.hartwig.hmftools.patientreporter.report.pages.AlterationEvidencePage;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.exception.DRException;

final class EvidenceReport {

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generate(@NotNull final SequencedPatientReport report) throws IOException, DRException {
        final EvidenceReportData evidenceReportData = ImmutableEvidenceReportData.of(report.civicAlterations());

        //MIVO: can't use multiPageList here because it does not pass its data source to child pages
        final ComponentBuilder<?, ?> reportPages = cmp.verticalList()
                .add(AlterationEvidencePage.reportComponent(report.sampleReport()))
                .add(cmp.pageBreak())
                .add(AlterationDebugPage.reportComponent());

        return report().addDetail(reportPages).setDataSource(evidenceReportData.toDataSource());
    }
}

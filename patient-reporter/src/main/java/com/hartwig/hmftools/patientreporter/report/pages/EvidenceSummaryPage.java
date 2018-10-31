package com.hartwig.hmftools.patientreporter.report.pages;

import static com.hartwig.hmftools.patientreporter.report.Commons.HEADER_TO_TABLE_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.linkStyle;
import static com.hartwig.hmftools.patientreporter.report.Commons.monospaceBaseTable;
import static com.hartwig.hmftools.patientreporter.report.Commons.sectionHeaderStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;

import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.report.Commons;
import com.hartwig.hmftools.patientreporter.report.components.MainPageTopSection;
import com.hartwig.hmftools.patientreporter.report.components.SummaryPart;
import com.hartwig.hmftools.patientreporter.report.data.ClinicalTrialDataSource;
import com.hartwig.hmftools.patientreporter.report.data.EvidenceItemDataSource;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;

@Value.Immutable
@Value.Style(passAnnotations = NotNull.class,
             allParameters = true)
public abstract class EvidenceSummaryPage {

    @NotNull
    abstract AnalysedPatientReport report();

    @NotNull
    public ComponentBuilder<?, ?> reportComponent() {
        return cmp.verticalList(MainPageTopSection.buildWithImpliedPurity(Commons.TITLE_SEQUENCE,
                report().sampleReport()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                SummaryPart.summaryData(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                evidenceItemReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                clinicalTrialReport(report()));
    }

    @NotNull
    private static ComponentBuilder<?, ?> evidenceItemReport(@NotNull AnalysedPatientReport report) {
        final ComponentBuilder<?, ?> table = report.clinicalEvidence().size() > 0
                ? cmp.subreport(monospaceBaseTable().fields(EvidenceItemDataSource.evidenceItemFields())
                .columns(col.column("Event", EvidenceItemDataSource.EVENT_FIELD).setFixedWidth(120),
                        col.column("Drug", EvidenceItemDataSource.DRUG_FIELD),
                        col.column("Level", EvidenceItemDataSource.LEVEL_FIELD),
                        col.column("Response", EvidenceItemDataSource.RESPONSE_FIELD),
                        col.column("Source", EvidenceItemDataSource.SOURCE_FIELD)
                                .setHyperLink(hyperLink(EvidenceItemDataSource.sourceHyperlink()))
                                .setStyle(linkStyle()))
                .setDataSource(EvidenceItemDataSource.fromEvidenceItems(report.clinicalEvidence())))
                : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Clinical Evidence").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_VERTICAL_GAP),
                table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> clinicalTrialReport(@NotNull AnalysedPatientReport report) {
        final ComponentBuilder<?, ?> table = report.clinicalTrials().size() > 0
                ? cmp.subreport(monospaceBaseTable().fields(ClinicalTrialDataSource.clinicalTrialFields())
                .columns(col.column("Event", ClinicalTrialDataSource.EVENT_FIELD).setFixedWidth(120),
                        col.column("Trial", ClinicalTrialDataSource.TRIAL_FIELD),
                        col.column("Source", ClinicalTrialDataSource.SOURCE_FIELD)
                                .setHyperLink(hyperLink(ClinicalTrialDataSource.sourceHyperlink()))
                                .setStyle(linkStyle()),
                        col.column("CCMO", ClinicalTrialDataSource.CCMO_FIELD))
                .setDataSource(ClinicalTrialDataSource.fromClinicalTrials(report.clinicalTrials())))
                : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Clinical Trials").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_VERTICAL_GAP),
                table);
    }
}
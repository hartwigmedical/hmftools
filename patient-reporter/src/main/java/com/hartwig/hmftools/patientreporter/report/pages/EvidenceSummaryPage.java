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

import java.util.List;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.patientreporter.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.report.Commons;
import com.hartwig.hmftools.patientreporter.report.components.GenomicSummarySection;
import com.hartwig.hmftools.patientreporter.report.components.MainPageTopSection;
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
        return cmp.verticalList(MainPageTopSection.build(Commons.TITLE_SEQUENCE, report().sampleReport()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                GenomicSummarySection.build(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                tumorTypeSpecificEvidence(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                clinicalTrialReport(report()),
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                offLabelClinicalEvidence(report()));
    }

    @NotNull
    private static ComponentBuilder<?, ?> tumorTypeSpecificEvidence(@NotNull AnalysedPatientReport report) {
        return evidenceItemReport(report.tumorSpecificEvidence(), "Tumor Type Specific Evidence");
    }

    @NotNull
    private static ComponentBuilder<?, ?> offLabelClinicalEvidence(@NotNull AnalysedPatientReport report) {
        return evidenceItemReport(report.offLabelEvidence(), "Evidence On Other Tumor Types");
    }

    @NotNull
    private static ComponentBuilder<?, ?> evidenceItemReport(@NotNull List<EvidenceItem> items, @NotNull String title) {
        final ComponentBuilder<?, ?> table = items.size() > 0
                ? cmp.subreport(monospaceBaseTable().fields(EvidenceItemDataSource.evidenceItemFields())
                .columns(col.column("Event", EvidenceItemDataSource.EVENT_FIELD).setMinWidth(140),
                        col.column("Match", EvidenceItemDataSource.SCOPE_FIELD),
                        col.column("Drug", EvidenceItemDataSource.DRUG_FIELD).setMinWidth(190),
                        col.column("Level", EvidenceItemDataSource.LEVEL_FIELD),
                        col.column("Response", EvidenceItemDataSource.RESPONSE_FIELD),
                        col.column("Source", EvidenceItemDataSource.SOURCE_FIELD)
                                .setHyperLink(hyperLink(EvidenceItemDataSource.sourceHyperlink()))
                                .setStyle(linkStyle()))
                .setDataSource(EvidenceItemDataSource.fromEvidenceItems(items)))
                : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text(title).setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_VERTICAL_GAP),
                table);
    }

    @NotNull
    private static ComponentBuilder<?, ?> clinicalTrialReport(@NotNull AnalysedPatientReport report) {
        final ComponentBuilder<?, ?> table = report.clinicalTrials().size() > 0
                ? cmp.subreport(monospaceBaseTable().fields(ClinicalTrialDataSource.clinicalTrialFields())
                .columns(col.column("Event", ClinicalTrialDataSource.EVENT_FIELD).setMinWidth(140),
                        col.column("Trial", ClinicalTrialDataSource.TRIAL_FIELD),
                        col.column("Source", ClinicalTrialDataSource.SOURCE_FIELD)
                                .setHyperLink(hyperLink(ClinicalTrialDataSource.sourceHyperlink()))
                                .setStyle(linkStyle()),
                        col.column("CCMO", ClinicalTrialDataSource.CCMO_FIELD))
                .setDataSource(ClinicalTrialDataSource.fromClinicalTrials(report.clinicalTrials())))
                : cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER));

        return cmp.verticalList(cmp.text("Clinical Trials (NL)").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(HEADER_TO_TABLE_VERTICAL_GAP),
                table);
    }
}
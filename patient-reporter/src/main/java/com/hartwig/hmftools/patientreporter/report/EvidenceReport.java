package com.hartwig.hmftools.patientreporter.report;

import static com.hartwig.hmftools.patientreporter.report.Commons.SECTION_VERTICAL_GAP;
import static com.hartwig.hmftools.patientreporter.report.Commons.fontStyle;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.exp;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;

import java.awt.Color;
import java.io.IOException;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.report.components.MainPageTopSection;
import com.hartwig.hmftools.patientreporter.report.data.Alteration;
import com.hartwig.hmftools.patientreporter.report.data.AlterationEvidence;
import com.hartwig.hmftools.patientreporter.report.data.AlterationMatch;
import com.hartwig.hmftools.patientreporter.report.data.EvidenceReportData;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.SubreportBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.exception.DRException;

class EvidenceReport {
    private static final int PADDING = 3;

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generate(@NotNull final PatientReport report, @NotNull final HmfReporterData reporterData)
            throws IOException, DRException {

        final EvidenceReportData evidenceReportData =
                EvidenceReportData.of(report, reporterData.geneModel(), reporterData.doidMapping().doidsForTumorType(report.tumorType()));

        // @formatter:off
        final VerticalListBuilder totalReport = cmp.verticalList(
                    MainPageTopSection.build("HMF Civic Evidence Supplement", report),
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    cmp.text("Knowledgebase drug association of reported genomic alterations").setStyle(sectionHeaderStyle().setFontSize(15)),
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    conciseEvidenceSection(),
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    cmp.text("Civic matching variants").setStyle(sectionHeaderStyle()),
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    civicMatchingVariantsTable(),
                    cmp.verticalGap(SECTION_VERTICAL_GAP));
        // @formatter:on

        return report().addDetail(totalReport).setDataSource(evidenceReportData.toDataSource());
    }

    @NotNull
    private static ComponentBuilder<?, ?> conciseEvidenceSection() throws IOException, DRException {
        return cmp.horizontalList(cmp.horizontalGap(20), alterationEvidenceTable(), cmp.horizontalGap(20));
    }

    @NotNull
    private static ComponentBuilder<?, ?> alterationEvidenceTable() throws IOException, DRException {
        final int ALTERATION_WIDTH = 150;
        final int SIGNIFICANCE_WIDTH = 100;
        final int DRUGS_WIDTH = 200;
        final int SOURCE_WIDTH = 70;

        //@formatter:off
        final SubreportBuilder subtable = cmp.subreport(
                baseTable().setColumnStyle(dataStyle()).fields(AlterationEvidence.SOURCE_URL)
                    .columns(
                        col.column(AlterationEvidence.SIGNIFICANCE).setFixedWidth(SIGNIFICANCE_WIDTH).setMinHeight(25),
                        col.column(AlterationEvidence.DRUGS).setFixedWidth(DRUGS_WIDTH),
                        col.column(AlterationEvidence.SOURCE).setHyperLink(hyperLink(AlterationEvidence.sourceHyperlink()))
                                .setStyle(dataLinkStyle()).setFixedWidth(SOURCE_WIDTH)))
                .setDataSource(exp.subDatasourceBeanCollection("evidence"));

        final ComponentBuilder<?, ?> tableHeader = cmp.horizontalList(
                cmp.text("Alteration").setStyle(tableHeaderStyle()).setFixedWidth(ALTERATION_WIDTH),
                cmp.text("Significance").setStyle(tableHeaderStyle()).setFixedWidth(SIGNIFICANCE_WIDTH),
                cmp.text("Association(Lv)").setStyle(tableHeaderStyle()).setFixedWidth(DRUGS_WIDTH),
                cmp.text("Source").setStyle(tableHeaderStyle()).setFixedWidth(SOURCE_WIDTH));

        return cmp.subreport(
                baseTable().setColumnStyle(dataStyle()).title(tableHeader)
                    .columns(
                        col.column(Alteration.ALTERATION).setFixedWidth(ALTERATION_WIDTH),
                        col.componentColumn(subtable))
                .noData(cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))))
                .setDataSource(exp.subDatasourceBeanCollection("alterationsWithEvidence"));
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> civicMatchingVariantsTable() throws IOException, DRException {
        final int ALTERATION_WIDTH = 75;

        //@formatter:off
        final SubreportBuilder subtable = cmp.subreport(
                baseTable().setColumnStyle(dataStyle()).fields(AlterationMatch.SUMMARY_URL)
                    .columns(
                        col.column(AlterationMatch.MATCH_TYPE).setFixedWidth(40).setMinHeight(25),
                        col.column(AlterationMatch.NAME).setHyperLink(hyperLink(AlterationMatch.civicSummaryHyperlink()))
                                .setStyle(dataLinkStyle()).setFixedWidth(80),
                        col.column(AlterationMatch.VARIANT_TYPE).setFixedWidth(80),
                        col.column(AlterationMatch.CHROMOSOME).setFixedWidth(30),
                        col.column(AlterationMatch.START).setFixedWidth(45),
                        col.column(AlterationMatch.STOP).setFixedWidth(45),
                        col.column(AlterationMatch.CHROMOSOME2).setFixedWidth(30),
                        col.column(AlterationMatch.START2).setFixedWidth(30),
                        col.column(AlterationMatch.STOP2).setFixedWidth(30),
                        col.column(AlterationMatch.HGVS_EXPRESSIONS)))
                .setDataSource(exp.subDatasourceBeanCollection("matches"));

        final ComponentBuilder<?, ?> tableHeader = cmp.horizontalList(
                cmp.text("Alteration").setStyle(tableHeaderStyle()).setFixedWidth(ALTERATION_WIDTH),
                cmp.text("Match").setStyle(tableHeaderStyle()).setFixedWidth(40),
                cmp.text("Name").setStyle(tableHeaderStyle()).setFixedWidth(80),
                cmp.text("Type").setStyle(tableHeaderStyle()).setFixedWidth(80),
                cmp.text("Chr").setStyle(tableHeaderStyle()).setFixedWidth(30),
                cmp.text("Start").setStyle(tableHeaderStyle()).setFixedWidth(45),
                cmp.text("Stop").setStyle(tableHeaderStyle()).setFixedWidth(45),
                cmp.text("Chr2").setStyle(tableHeaderStyle()).setFixedWidth(30),
                cmp.text("Start2").setStyle(tableHeaderStyle()).setFixedWidth(30),
                cmp.text("Stop2").setStyle(tableHeaderStyle()).setFixedWidth(30),
                cmp.text("Hgvs Expressions").setStyle(tableHeaderStyle()));

        return cmp.subreport(
                baseTable().setColumnStyle(dataStyle()).title(tableHeader)
                    .columns(
                        col.column(Alteration.ALTERATION).setFixedWidth(ALTERATION_WIDTH),
                        col.componentColumn(subtable))
                .noData(cmp.text("None").setStyle(fontStyle().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))))
                .setDataSource(exp.subDatasourceBeanCollection("alterations"));
        // @formatter:on
    }

    @NotNull
    private static JasperReportBuilder baseTable() {
        return report().setColumnStyle(dataStyle()).setColumnTitleStyle(tableHeaderStyle());
    }

    @NotNull
    private static StyleBuilder sectionHeaderStyle() {
        return fontStyle().bold().setFontSize(12).setHorizontalTextAlignment(HorizontalTextAlignment.CENTER);
    }

    @NotNull
    private static StyleBuilder tableHeaderStyle() {
        return Commons.tableHeaderStyle().setFontSize(9).setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder dataStyle() {
        return Commons.smallDataTableStyle().setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder dataLinkStyle() {
        return dataStyle().setForegroundColor(Color.BLUE);
    }
}

package com.hartwig.hmftools.patientreporter.report.layout;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.exp;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.patientreporter.report.data.AlterationEvidenceReporterData;
import com.hartwig.hmftools.patientreporter.report.data.AlterationMatch;
import com.hartwig.hmftools.patientreporter.report.data.AlterationReporterData;
import com.hartwig.hmftools.patientreporter.report.data.VariantReporterData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.SubreportBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.exception.DRException;

public class EvidenceLayout {

    private static final Logger LOGGER = LogManager.getLogger(EvidenceLayout.class);

    private static final String FONT = "Times New Roman";
    private static final Color BORKIE_COLOR = new Color(221, 235, 247);
    private static final int SECTION_VERTICAL_GAP = 25;
    private static final int PADDING = 3;

    @NotNull
    private final String reportDirectory;

    public EvidenceLayout(@NotNull final String reportDirectory) {
        this.reportDirectory = reportDirectory;
    }

    //    @Override
    public void writeSequenceReport(@NotNull final PatientReport report, @NotNull final HmfReporterData reporterData)
            throws IOException, DRException {
        final JasperReportBuilder reportBuilder = generatePatientReport(report, reporterData);
        writeReport(report.sample(), reportBuilder);
    }

    private void writeReport(@NotNull final String sample, @NotNull final JasperReportBuilder report)
            throws FileNotFoundException, DRException {
        final String fileName = fileName(sample);
        if (Files.exists(new File(fileName).toPath())) {
            LOGGER.warn(" Could not write report as it already exists: " + fileName);
        } else {
            report.toPdf(new FileOutputStream(fileName));
            LOGGER.info(" Created evidence report at " + fileName);
        }
    }

    @NotNull
    private String fileName(@NotNull final String sample) {
        return reportDirectory + File.separator + sample + "_evidence_report.pdf";
    }

    @VisibleForTesting
    @NotNull
    public static JasperReportBuilder generatePatientReport(@NotNull final PatientReport report,
            @NotNull final HmfReporterData reporterData) throws IOException, DRException {

        final VariantReporterData variantReporterData =
                VariantReporterData.of(report, reporterData.geneModel(), reporterData.doidMapping().doidsForTumorType(report.tumorType()));

        // @formatter:off
        final VerticalListBuilder totalReport = cmp.verticalList(
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    cmp.text("Knowledgebase drug association of reported genomic alterations").setStyle(sectionHeaderStyle().setFontSize(15)),
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    conciseEvidenceSection(),
                    cmp.verticalGap(SECTION_VERTICAL_GAP),
                    cmp.text("Civic matching variants"),
                    civicMatchingVariantsTable(),
                    cmp.verticalGap(SECTION_VERTICAL_GAP));

        final VerticalListBuilder noDataReport = cmp.verticalList()
                .add(cmp.verticalGap(SECTION_VERTICAL_GAP),
                    cmp.text("HMF Sequencing Report v" + PatientReporterApplication.VERSION + " - Could not find any civic evidence items").setStyle(sectionHeaderStyle()),
                    cmp.verticalGap(SECTION_VERTICAL_GAP));
        // @formatter:on

        return report().addDetail(totalReport).noData(noDataReport).setDataSource(variantReporterData.toDataSource());
    }

    @NotNull
    private static ComponentBuilder<?, ?> conciseEvidenceSection() throws IOException, DRException {
        return cmp.horizontalList(cmp.horizontalGap(20), alterationEvidenceTable(), cmp.horizontalGap(20));
    }

    @NotNull
    private static ComponentBuilder<?, ?> alterationEvidenceTable() throws IOException, DRException {
        final int ALTERATION_WIDTH = 135;
        final int SIGNIFICANCE_WIDTH = 70;
        final int DRUGS_WIDTH = 170;
        final int SOURCE_WIDTH = 60;

        //@formatter:off
        final SubreportBuilder subtable = cmp.subreport(
                baseTable().setColumnStyle(dataStyle())
                    .columns(
                        col.column(AlterationEvidenceReporterData.SIGNIFICANCE).setFixedWidth(SIGNIFICANCE_WIDTH).setMinHeight(25),
                        col.column(AlterationEvidenceReporterData.DRUGS).setFixedWidth(DRUGS_WIDTH),
                        col.column(AlterationEvidenceReporterData.SOURCE).setFixedWidth(SOURCE_WIDTH)))
                .setDataSource(exp.subDatasourceBeanCollection("evidence"));

        final ComponentBuilder<?, ?> tableHeader = cmp.horizontalList(
                cmp.text("Alteration").setStyle(tableHeaderStyle()).setFixedWidth(ALTERATION_WIDTH),
                cmp.text("Significance").setStyle(tableHeaderStyle()).setFixedWidth(SIGNIFICANCE_WIDTH),
                cmp.text("Association(Lv)").setStyle(tableHeaderStyle()).setFixedWidth(DRUGS_WIDTH),
                cmp.text("Source").setStyle(tableHeaderStyle()).setFixedWidth(SOURCE_WIDTH));

        return cmp.subreport(
                baseTable().setColumnStyle(dataStyle()).title(tableHeader)
                    .columns(
                        col.column(AlterationReporterData.ALTERATION).setFixedWidth(ALTERATION_WIDTH),
                        col.componentColumn(subtable)))
                .setDataSource(exp.subDatasourceBeanCollection("alterations"));
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
                                .setStyle(linkStyle()).setFixedWidth(80),
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
                        col.column(AlterationReporterData.ALTERATION).setFixedWidth(ALTERATION_WIDTH),
                        col.componentColumn(subtable)))
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
        return fontStyle().bold()
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setFontSize(9)
                .setBorder(stl.pen1Point())
                .setBackgroundColor(BORKIE_COLOR)
                .setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder dataStyle() {
        return fontStyle().setFontSize(7)
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setBorder(stl.penThin().setLineColor(Color.black))
                .setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder linkStyle() {
        return dataStyle().setForegroundColor(Color.BLUE);
    }

    @NotNull
    private static StyleBuilder fontStyle() {
        return stl.style().setFontName(FONT);
    }
}

package com.hartwig.hmftools.patientreporter.report.layout;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.exp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.patientreporter.report.data.CivicVariantReporterData;
import com.hartwig.hmftools.patientreporter.report.data.EvidenceItemReporterData;
import com.hartwig.hmftools.patientreporter.report.data.VariantReporterData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.MultiPageListBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.exception.DRException;
import net.sf.jasperreports.engine.data.JRBeanCollectionDataSource;

public class EvidenceLayout {

    private static final Logger LOGGER = LogManager.getLogger(EvidenceLayout.class);

    private static final String FONT = "Times New Roman";
    private static final Color BORKIE_COLOR = new Color(221, 235, 247);

    private static final int TEXT_HEADER_INDENT = 30;
    private static final int TEXT_DETAIL_INDENT = 40;
    private static final int LIST_INDENT = 5;
    private static final int HEADER_TO_DETAIL_VERTICAL_GAP = 8;
    private static final int DETAIL_TO_DETAIL_VERTICAL_GAP = 4;
    private static final int SECTION_VERTICAL_GAP = 25;
    private static final int PADDING = 1;

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
            LOGGER.info(" Created patient report at " + fileName);
        }
    }

    @NotNull
    private String fileName(@NotNull final String sample) {
        return reportDirectory + File.separator + sample + "_evidence_report.pdf";
    }

    @VisibleForTesting
    @NotNull
    public static JasperReportBuilder generatePatientReport(@NotNull final PatientReport report,
            @NotNull final HmfReporterData reporterData) {
        // @formatter:off

        final MultiPageListBuilder totalReport = cmp.multiPageList().add(
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("HMF Sequencing Report v" + PatientReporterApplication.VERSION + " - Civic Evidence Items").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP));

        final ComponentBuilder<?, ?> evidenceItemsPage = cmp.verticalList(evidenceSection());
        // @formatter:on

        final List<VariantReporterData> variantReporterData = VariantReporterData.of(report, reporterData);

        return report().addDetail(evidenceItemsPage).setDataSource(new JRBeanCollectionDataSource(Lists.newArrayList(variantReporterData)));
    }

    @NotNull
    private static ComponentBuilder<?, ?> evidenceSection() {
        return cmp.verticalList(cmp.text(VariantReporterData.VARIANT_NAME),
                cmp.subreport(report().detail(evidenceTable())).setDataSource(exp.subDatasourceBeanCollection("variants")));
    }

    @NotNull
    private static ComponentBuilder<?, ?> evidenceTable() {
        //@formatter:off
        final int fontSize = 7;
        return cmp.verticalList(
                cmp.text(CivicVariantReporterData.VARIANT).setStyle(tableHeaderStyle()),
                cmp.subreport(
                    baseTable().setColumnStyle(dataStyle().setFontSize(fontSize))
                        .columns(
                            col.column("Tumor Type", EvidenceItemReporterData.TUMOR_TYPE_FIELD).setFixedWidth(100),
                            col.column("Level", EvidenceItemReporterData.LEVEL_FIELD).setFixedWidth(50),
                            col.column("Direction", EvidenceItemReporterData.DIRECTION_FIELD).setFixedWidth(100),
                            col.column("Significance", EvidenceItemReporterData.SIGNIFICANCE_FIELD).setFixedWidth(100),
                            col.column("Drugs", EvidenceItemReporterData.DRUGS_FIELD)))
                    .setDataSource(exp.subDatasourceBeanCollection("evidenceItems")),
                cmp.verticalGap(10));
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> conciseEvidenceSection() {
        return cmp.verticalList(cmp.text(CivicVariantReporterData.VARIANT), conciseEvidenceTable());
    }

    @NotNull
    private static ComponentBuilder<?, ?> conciseEvidenceTable() {
        //@formatter:off
        final int fontSize = 7;
        return cmp.subreport(
                baseTable().setColumnStyle(dataStyle().setFontSize(fontSize))
                    .title(cmp.text(CivicVariantReporterData.VARIANT))
                    .columns(
                        col.column("Significance", EvidenceItemReporterData.SIGNIFICANCE_FIELD).setFixedWidth(100),
                        col.column("Drugs", EvidenceItemReporterData.DRUGS_FIELD).setFixedWidth(400)))
                .setDataSource(exp.subDatasourceBeanCollection("evidenceItems"));
        // @formatter:on
    }

    @NotNull
    private static JasperReportBuilder baseTable() {
        return report().setColumnStyle(dataStyle()).setColumnTitleStyle(tableHeaderStyle()).highlightDetailEvenRows();
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
                .setFontSize(10)
                .setBorder(stl.pen1Point())
                .setBackgroundColor(BORKIE_COLOR)
                .setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder dataStyle() {
        return fontStyle().setFontSize(8)
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
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

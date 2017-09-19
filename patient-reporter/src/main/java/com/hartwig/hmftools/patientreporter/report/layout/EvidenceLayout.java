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
import com.hartwig.hmftools.patientreporter.report.data.AlterationEvidenceReporterData;
import com.hartwig.hmftools.patientreporter.report.data.AlterationReporterData;
import com.hartwig.hmftools.patientreporter.report.data.CivicVariantReporterData;
import com.hartwig.hmftools.patientreporter.report.data.EvidenceItemReporterData;
import com.hartwig.hmftools.patientreporter.report.data.ImmutableAlterationEvidenceReporterData;
import com.hartwig.hmftools.patientreporter.report.data.ImmutableAlterationReporterData;
import com.hartwig.hmftools.patientreporter.report.data.VariantReporterData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.MultiPageListBuilder;
import net.sf.dynamicreports.report.builder.component.SubreportBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.exception.DRException;
import net.sf.jasperreports.engine.data.JRBeanCollectionDataSource;

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
            @NotNull final HmfReporterData reporterData) throws IOException, DRException {

        final List<VariantReporterData> variantReporterData = VariantReporterData.of(report, reporterData.geneModel());

        //TODO: get actual data
        final List<Object> conciseReporterData = Lists.newArrayList(new Object());
        final ComponentBuilder<?, ?> alterationsReport = cmp.subreport(report().addDetail(conciseEvidenceSection())
                .setDataSource(new JRBeanCollectionDataSource(Lists.newArrayList(conciseReporterData))));

        final ComponentBuilder<?, ?> evidenceItemsReport = cmp.subreport(report().addDetail(evidenceSection())
                .setDataSource(new JRBeanCollectionDataSource(Lists.newArrayList(variantReporterData))));

        // @formatter:off
        final MultiPageListBuilder totalReport = cmp.multiPageList().add(
                cmp.verticalGap(SECTION_VERTICAL_GAP),
                cmp.text("HMF Sequencing Report v" + PatientReporterApplication.VERSION + " - Civic Evidence Items").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(SECTION_VERTICAL_GAP))
                .add(alterationsReport)
                .newPage()
                .add(evidenceItemsReport);
        // @formatter:on

        return report().addDetail(totalReport).setDataSource(new JRBeanCollectionDataSource(Lists.newArrayList(variantReporterData)));
    }

    @NotNull
    private static ComponentBuilder<?, ?> evidenceSection() {
        return cmp.horizontalList(cmp.horizontalGap(20), cmp.verticalList(cmp.text(VariantReporterData.VARIANT_NAME),
                cmp.subreport(report().detail(evidenceTable())).setDataSource(exp.subDatasourceBeanCollection("variants"))),
                cmp.horizontalGap(20));
    }

    @NotNull
    private static ComponentBuilder<?, ?> evidenceTable() {
        //@formatter:off
        return cmp.verticalList(
                cmp.text(CivicVariantReporterData.VARIANT).setStyle(tableHeaderStyle()),
                cmp.subreport(
                    baseTable().setColumnStyle(dataStyle())
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
    private static ComponentBuilder<?, ?> conciseEvidenceSection() throws IOException, DRException {
        return cmp.horizontalList(cmp.horizontalGap(20), alterationEvidenceTable(), cmp.horizontalGap(20));
    }

    @NotNull
    private static ComponentBuilder<?, ?> alterationEvidenceTable() throws IOException, DRException {
        //@formatter:off
        final List<AlterationReporterData> alterationReporterData = Lists.newArrayList(
                ImmutableAlterationReporterData.of("EGFR", "p.Glu746_Ala750del", "", "", Lists.newArrayList(ImmutableAlterationEvidenceReporterData.of("sensitive", "erlotinib(A), afatinib(B), gefitinib(A), dacomitinib(B)", "CIViC"))),
                ImmutableAlterationReporterData.of("EGFR", "p.Thr790Met", "", "yes", Lists.newArrayList(
                        ImmutableAlterationEvidenceReporterData.of("resistant", "erlotinib(B), afatinib(B), gefitinib(B), dacomitinib(B)", "CIViC"),
                        ImmutableAlterationEvidenceReporterData.of("sensitive", "osimertinib(A)", "CIViC"))));

        final SubreportBuilder significanceSubreport = cmp.subreport(report().setColumnStyle(dataStyle()).columns(col.column(AlterationEvidenceReporterData.SIGNIFICANCE).setMinHeight(25))).setDataSource(exp.subDatasourceBeanCollection("evidence"));
        final SubreportBuilder drugsSubreport = cmp.subreport(report().setColumnStyle(dataStyle()).columns(col.column(AlterationEvidenceReporterData.DRUGS).setMinHeight(25))).setDataSource(exp.subDatasourceBeanCollection("evidence"));
        final SubreportBuilder sourceSubreport = cmp.subreport(report().setColumnStyle(dataStyle()).columns(col.column(AlterationEvidenceReporterData.SOURCE).setMinHeight(25))).setDataSource(exp.subDatasourceBeanCollection("evidence"));

        return cmp.subreport(
                baseTable().setColumnStyle(dataStyle())
                    .columns(
                        col.column("Alteration", AlterationReporterData.ALTERATION).setFixedWidth(135),
                        col.componentColumn("Significance", significanceSubreport).setStyle(dataStyle()).setFixedWidth(70),
                        col.componentColumn("Association(Lv)", drugsSubreport).setFixedWidth(170),
                        col.componentColumn("Source", sourceSubreport).setFixedWidth(60),
                        col.column("LOH", AlterationReporterData.LOH).setFixedWidth(40),
                        col.column("Subclonal", AlterationReporterData.SUBCLONAL).setFixedWidth(60)))
                .setDataSource(new JRBeanCollectionDataSource(Lists.newArrayList(alterationReporterData)));
//                .setDataSource(exp.subDatasourceBeanCollection("conciseEvidenceItems"));
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
                .setBorder(stl.penThin())
                .setPadding(PADDING);
    }

    @NotNull
    private static StyleBuilder fontStyle() {
        return stl.style().setFontName(FONT);
    }
}

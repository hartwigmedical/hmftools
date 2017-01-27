package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.patientreporter.PatientReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.column.TextColumnBuilder;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriter {

    private static final String DISCLAIMER = "This test is not performed in a diagnostic laboratory and "
            + "can therefore not be used to directly affect clinical decision making.";
    private static final String REF_TO_EXPLANATIONS =
            "Short explanations on the various fields can be " + "found on the final page of this report.";

    @NotNull
    private final String outputDirectory;
    @NotNull
    private final String hmfLogo;

    public PDFWriter(@NotNull final String outputDirectory, @NotNull final String hmfLogo) {
        this.outputDirectory = outputDirectory;
        this.hmfLogo = hmfLogo;
    }

    @NotNull
    public String write(@NotNull final PatientReport report) throws FileNotFoundException, DRException {
        final String fileName = outputDirectory + File.separator + report.sample() + "_hmf_report.pdf";
        final JasperReportBuilder jasperReportBuilder = generatePatientReport(report, hmfLogo);

        jasperReportBuilder.toPdf(new FileOutputStream(fileName));

        return fileName;
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generatePatientReport(@NotNull final PatientReport report,
            @NotNull final String hmfLogoPath) {

        // @formatter:off
        final ComponentBuilder<?,?> disclaimer = cmp.verticalList(
                cmp.horizontalList(
                        cmp.horizontalGap(30),
                        cmp.text("About this report").setStyle(fontStyle().bold().setFontSize(11))),
                cmp.verticalGap(5),
                cmp.horizontalList(
                        cmp.horizontalGap(30),
                        cmp.text("   1. " + DISCLAIMER).setStyle(fontStyle())),
                cmp.horizontalList(
                        cmp.horizontalGap(30),
                        cmp.text("   2. " + REF_TO_EXPLANATIONS).setStyle(fontStyle()))
        );

        return report().title(
                cmp.verticalList(
                        topSection(report, hmfLogoPath),
                        cmp.verticalGap(25),
                        disclaimer,
                        cmp.verticalGap(25),
                        variantReport(report),
                        cmp.verticalGap(25),
                        copyNumberReport(report))
                );
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> topSection(@NotNull final PatientReport report,
            @NotNull final String hmfLogoPath) {
        // @formatter:off
        final ComponentBuilder<?, ?> mainDiagnosisInfo = cmp.horizontalList(
                cmp.verticalList(
                        cmp.text("Report Date").setStyle(tableHeaderStyle()),
                        cmp.currentDate().setPattern("dd-MMM-yyyy").setStyle(dataTableStyle())),
                cmp.verticalList(
                        cmp.text("Tumor Type").setStyle(tableHeaderStyle()),
                        cmp.text(report.tumorType()).setStyle(dataTableStyle()))
        );

        return cmp.horizontalList(
                cmp.image(hmfLogoPath),
                cmp.verticalList(
                        cmp.text("HMF Sequencing Report - " + report.sample()).
                                setStyle(fontStyle().bold().setFontSize(14)
                                        .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE))
                                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                                .setHeight(50),
                        mainDiagnosisInfo)
        );
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> variantReport(@NotNull final PatientReport report) {
        // @formatter:off
        return cmp.verticalList(
                cmp.text("Somatic Variants").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(6),
                cmp.subreport(baseTable().fields(PatientDataSource.variantFields())
                        .columns(
                            col.column("Gene", PatientDataSource.GENE_FIELD),
                            col.column("Position", PatientDataSource.POSITION_FIELD),
                            col.column("Variant", PatientDataSource.VARIANT_FIELD),
                            transcriptColumn(),
                            col.componentColumn("Effect", effectColumn()),
                            col.column("Cosmic", PatientDataSource.COSMIC_FIELD)
                                    .setHyperLink(hyperLink(new COSMICLinkExpression())).setStyle(linkStyle()),
                            col.column("VAF", PatientDataSource.ALLELE_FREQUENCY_FIELD)))
                        .setDataSource(PatientDataSource.fromVariants(report.variants())),
                cmp.verticalGap(15),
                cmp.text("Mutational Load: " + Integer.toString(report.mutationalLoad())).setStyle(tableHeaderStyle())
        );
        // @formatter:on
    }

    @NotNull
    private static ComponentBuilder<?, ?> copyNumberReport(@NotNull final PatientReport report) {
        // @formatter:off
        return cmp.verticalList(
                cmp.text("Somatic Copy Numbers").setStyle(sectionHeaderStyle()),
                cmp.verticalGap(6),
                cmp.subreport(baseTable().fields(PatientDataSource.copyNumberFields())
                        .columns(
                            col.column("Gene", PatientDataSource.GENE_FIELD),
                            transcriptColumn(),
                            col.column("Copies", PatientDataSource.COPY_NUMBER_FIELD))
                        .setDataSource(PatientDataSource.fromCopyNumbers(report.copyNumbers())))
        );
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
        return fontStyle().bold().setHorizontalTextAlignment(HorizontalTextAlignment.CENTER).setFontSize(10).setBorder(
                stl.pen1Point()).setBackgroundColor(new Color(210, 210, 210));
    }

    @NotNull
    private static StyleBuilder dataStyle() {
        return fontStyle().setFontSize(8).setHorizontalTextAlignment(
                HorizontalTextAlignment.CENTER).setVerticalTextAlignment(VerticalTextAlignment.MIDDLE);
    }

    @NotNull
    private static StyleBuilder dataTableStyle() {
        return dataStyle().setBorder(stl.pen1Point());
    }

    @NotNull
    private static StyleBuilder linkStyle() {
        return dataStyle().setForegroundColor(Color.BLUE);
    }

    @NotNull
    private static StyleBuilder fontStyle() {
        return stl.style().setFontName("Times New Roman");
    }

    @NotNull
    private static TextColumnBuilder<?> transcriptColumn() {
        return col.column("Transcript", PatientDataSource.TRANSCRIPT_FIELD).setWidth(150).setHyperLink(
                hyperLink(new TranscriptLinkExpression())).setStyle(linkStyle());
    }

    @NotNull
    private static ComponentBuilder<?, ?> effectColumn() {
        return cmp.verticalList(
                cmp.horizontalList(cmp.text(DataExpression.fromField(PatientDataSource.HGVS_CODING_FIELD)),
                        cmp.text(DataExpression.fromField(PatientDataSource.HGVS_PROTEIN_FIELD))),
                cmp.text(DataExpression.fromField(PatientDataSource.EFFECT_FIELD)));
    }

    private static class COSMICLinkExpression extends AbstractSimpleExpression<String> {
        public String evaluate(@NotNull final ReportParameters data) {
            return "http://grch37-cancer.sanger.ac.uk/cosmic/mutation/overview?id=" + data.getValue(
                    PatientDataSource.COSMIC_NR_FIELD.getName());
        }
    }

    private static class TranscriptLinkExpression extends AbstractSimpleExpression<String> {
        public String evaluate(@NotNull final ReportParameters data) {
            return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=" + data.getValue(
                    PatientDataSource.TRANSCRIPT_FIELD.getName());
        }
    }

}

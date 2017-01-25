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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.HorizontalListBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriter {

    private static final Logger LOGGER = LogManager.getLogger(PDFWriter.class);

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
        final String fileName = outputDirectory + File.separator + report.sample() + "_report.pdf";
        final JasperReportBuilder jasperReportBuilder = generatePatientReport(report, hmfLogo);

        // KODU: This line causes net.sf.jasperreports.engine.JRRuntimeException in production...
        jasperReportBuilder.toPdf(new FileOutputStream(fileName));

        return fileName;
    }

    @VisibleForTesting
    @NotNull
    static JasperReportBuilder generatePatientReport(@NotNull final PatientReport report,
            @NotNull final String hmfLogoPath) {
        // @formatter:off
        final StyleBuilder columnStyle = stl.style()
                .setFontSize(8)
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER);

        final StyleBuilder columnTitleStyle = stl.style()
                        .bold()
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                        .setBorder(stl.pen1Point())
                        .setBackgroundColor(Color.LIGHT_GRAY);
        final StyleBuilder linkStyle = stl.style(columnStyle)
                .setForegroundColor(Color.BLUE);

        final HorizontalListBuilder mainTitle =
                cmp.horizontalList().add(
                        cmp.image(hmfLogoPath),
                        cmp.text("HMF Sequencing Report - " + report.sample()).setStyle(stl.style()
                                .bold()
                                .setFontSize(12)
                                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)));

        final ComponentBuilder<?,?> variantAnnotationField = cmp.verticalList(
                cmp.horizontalList(
                        cmp.text(DataExpression.fromField(PatientDataSource.HGVS_CODING_FIELD)),
                        cmp.text(DataExpression.fromField(PatientDataSource.HGVS_PROTEIN_FIELD))),
                cmp.text(DataExpression.fromField(PatientDataSource.EFFECT_FIELD)));

        final JasperReportBuilder variantReport = report()
                .fields(PatientDataSource.variantFields())
                .columns(col.column("Gene", PatientDataSource.GENE_FIELD),
                         col.column("Position", PatientDataSource.POSITION_FIELD),
                         col.column("Variant", PatientDataSource.VARIANT_FIELD),
                         col.column("Transcript", PatientDataSource.TRANSCRIPT_FIELD).setWidth(150)
                                 .setHyperLink(hyperLink(new TranscriptLinkExpression()))
                                 .setStyle(linkStyle),
                         col.componentColumn("Annotation", variantAnnotationField),
                         col.column("Cosmic", PatientDataSource.COSMIC_FIELD)
                                 .setHyperLink(hyperLink(new COSMICLinkExpression()))
                                 .setStyle(linkStyle),
                         col.column("VAF", PatientDataSource.ALLELE_FREQUENCY_FIELD))
                .setColumnStyle(columnStyle)
                .setColumnTitleStyle(columnTitleStyle)
                .highlightDetailEvenRows();

        final JasperReportBuilder copyNumberReport = report()
                .fields(PatientDataSource.copyNumberFields())
                .columns(col.column("Gene", PatientDataSource.GENE_FIELD),
                         col.column("Transcript", PatientDataSource.TRANSCRIPT_FIELD).setWidth(150)
                                 .setHyperLink(hyperLink(new TranscriptLinkExpression()))
                                 .setStyle(linkStyle),
                         col.column("Number of copies?", PatientDataSource.COPY_NUMBER_FIELD))
                .setColumnStyle(columnStyle)
                .setColumnTitleStyle(columnTitleStyle)
                .highlightDetailEvenRows();

        return report().title(
                cmp.verticalList(
                        mainTitle,
                        cmp.verticalGap(20),
                        cmp.text("Mutational Load: " + Integer.toString(report.mutationalLoad()))
                            .setStyle(columnTitleStyle),
                        cmp.verticalGap(20),
                        cmp.subreport(variantReport)
                                .setDataSource(PatientDataSource.fromVariants(report.variants())),
                        cmp.verticalGap(20),
                        cmp.subreport(copyNumberReport)
                                .setDataSource(PatientDataSource.fromCopyNumbers(report.copyNumbers()))));
        // @formatter:on
    }

    private static class COSMICLinkExpression extends AbstractSimpleExpression<String> {
        public String evaluate(@NotNull final ReportParameters data) {
            return "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=" + data.getValue(
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

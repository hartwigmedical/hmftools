package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.hyperLink;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.component.ComponentBuilder;
import net.sf.dynamicreports.report.builder.component.HorizontalListBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final String HMF_LOGO = Resources.getResource("pdf/hartwig_logo.jpg").getPath();

    @Test
    @Ignore
    public void experimentWithReport() throws DRException, FileNotFoundException {
        final String sample = "CPCT11111111T";
        final List<VariantReport> variants = Lists.newArrayList(
                new VariantReport.Builder().gene("PIK3CA").position("3:178936082").ref("G").alt("A").transcript(
                        "ENST00000263967.3").hgvsCoding("c.1624G>A").hgvsProtein("p.Glu542Lys").consequence(
                        "missense variant").cosmicID("COSM125369").alleleReadCount(75).totalReadCount(151).build(),
                new VariantReport.Builder().gene("PIK3CA").position("3:178936082").ref("G").alt("A").transcript(
                        "ENST00000263967.3").hgvsCoding("c.1624G>A").hgvsProtein("p.Glu542Lys").consequence(
                        "missense variant").cosmicID("").alleleReadCount(75).totalReadCount(151).build(),
                new VariantReport.Builder().gene("PIK3CA").position("3:178936082").ref("G").alt("A").transcript(
                        "ENST00000263967.3").hgvsCoding("c.1624G>A").hgvsProtein("p.Glu542Lys").consequence(
                        "missense variant").cosmicID("COSM125369").alleleReadCount(75).totalReadCount(151).build());

        final List<CopyNumberReport> copyNumbers = Lists.newArrayList();
        final int mutationalLoad = 0;

        // @formatter:off
        final ComponentBuilder<?,?> annotationField = cmp.verticalList(
                cmp.horizontalList(
                        cmp.text(DataExpression.fromField(PatientDataSource.HGVS_CODING_FIELD)),
                        cmp.text(DataExpression.fromField(PatientDataSource.HGVS_PROTEIN_FIELD))),
                cmp.text(DataExpression.fromField(PatientDataSource.EFFECT_FIELD)));

        final StyleBuilder columnStyle = stl.style()
                .setFontSize(8)
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER);

        final StyleBuilder linkStyle = stl.style(columnStyle)
                .setForegroundColor(Color.BLUE);

        final OutputStream output = new FileOutputStream("/Users/kduyvesteyn/hmf/tmp/report.pdf");


        HorizontalListBuilder title =
                cmp.horizontalList().add(
                        cmp.image(HMF_LOGO),
                        cmp.text("HMF Sequencing Report - " + sample).setStyle(stl.style()
                                .bold()
                                .setFontSize(12)
                                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)));

        JasperReportBuilder variantTable = report()
                .fields(PatientDataSource.variantFields())
                .columns(col.column("Gene", PatientDataSource.GENE_FIELD),
                         col.column("Position", PatientDataSource.POSITION_FIELD),
                         col.column("Variant", PatientDataSource.VARIANT_FIELD),
                         col.column("Transcript", PatientDataSource.TRANSCRIPT_FIELD).setWidth(150)
                                 .setHyperLink(hyperLink(new TranscriptLinkExpression()))
                                 .setStyle(linkStyle),
                         col.componentColumn("Annotation", annotationField),
                         col.column("Cosmic", PatientDataSource.COSMIC_FIELD)
                                 .setHyperLink(hyperLink(new COSMICLinkExpression()))
                                 .setStyle(linkStyle),
                         col.column("VAF", PatientDataSource.ALLELE_FREQUENCY_FIELD))
                .setColumnStyle(columnStyle)
                .setColumnTitleStyle(stl.style()
                        .bold()
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                        .setBorder(stl.pen1Point())
                        .setBackgroundColor(Color.LIGHT_GRAY))
                .highlightDetailEvenRows()
                .setDataSource(PatientDataSource.fromVariants(variants));

        report().title(
                cmp.verticalList(
                        title,
                        cmp.subreport(variantTable).setDataSource(PatientDataSource.fromVariants(variants)),
                        cmp.verticalGap(20),
                        cmp.subreport(variantTable).setDataSource(PatientDataSource.fromVariants(variants))))
                .show()
                .print();
        // @formatter:on
    }

    private class COSMICLinkExpression extends AbstractSimpleExpression<String> {
        public String evaluate(@NotNull final ReportParameters data) {
            return "http://cancer.sanger.ac.uk/cosmic/mutation/overview?id=" + data.getValue(
                    PatientDataSource.COSMIC_NR_FIELD.getName());
        }
    }

    private class TranscriptLinkExpression extends AbstractSimpleExpression<String> {
        public String evaluate(@NotNull final ReportParameters data) {
            return "http://grch37.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=" + data.getValue(
                    PatientDataSource.TRANSCRIPT_FIELD.getName());
        }
    }
}
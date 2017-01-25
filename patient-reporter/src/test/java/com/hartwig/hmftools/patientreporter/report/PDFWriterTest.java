package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;
import static net.sf.dynamicreports.report.builder.DynamicReports.type;

import java.awt.Color;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.junit.Ignore;
import org.junit.Test;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.definition.ReportParameters;
import net.sf.dynamicreports.report.exception.DRException;
import net.sf.jasperreports.engine.JRDataSource;

public class PDFWriterTest {

    private static final String HMF_LOGO = Resources.getResource("pdf/hartwig_logo.png").getPath();

    @Test
    @Ignore
    public void experimentWithReport() throws DRException {
        final String sample = "CPCT11111111T";
        final List<VariantReport> variants = Lists.newArrayList(
                new VariantReport.Builder().gene("PIK3CA").position("3:178936082").ref("G").alt("A").transcript(
                        "ENST00000263967.3").hgvsCoding("c.1624G>A").hgvsProtein("p.Glu542Lys").consequence(
                        "missense variant").cosmicID("COSM125369").alleleReadCount(75).totalReadCount(151).build(),
                new VariantReport.Builder().gene("PIK3CA").position("3:178936082").ref("G").alt("A").transcript(
                        "ENST00000263967.3").hgvsCoding("c.1624G>A").hgvsProtein("p.Glu542Lys").consequence(
                        "missense variant").cosmicID("COSM125369").alleleReadCount(75).totalReadCount(151).build(),
                new VariantReport.Builder().gene("PIK3CA").position("3:178936082").ref("G").alt("A").transcript(
                        "ENST00000263967.3").hgvsCoding("c.1624G>A").hgvsProtein("p.Glu542Lys").consequence(
                        "missense variant").cosmicID("COSM125369").alleleReadCount(75).totalReadCount(151).build());

        final JRDataSource variantDataSource = DataSourceCreator.fromVariants(variants);

        final List<CopyNumberReport> copyNumbers = Lists.newArrayList();
        //        final int mutationalLoad = 0;
        //        final PatientReport patientReport = new PatientReport(sample, variants, copyNumbers, mutationalLoad);

        VerticalListBuilder listBuilder = cmp.verticalList(cmp.text(new Item1Expression()),
                cmp.text(new Item2Expression()));

        // @formatter:off
        report()
                .title(
                    cmp.horizontalList().add(
                        cmp.image(HMF_LOGO),
                        cmp.text("HMF Sequencing Report - " + sample).setStyle(stl.style()
                                .bold()
                                .setFontSize(12)
                                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE))))
                .fields(DataSourceCreator.variantFields())
                .columns(col.column("Gene", DataSourceCreator.GENE_FIELD, type.stringType()),
                         col.column("Position", DataSourceCreator.POSITION_FIELD, type.stringType()),
                         col.column("Variant", DataSourceCreator.VARIANT_FIELD, type.stringType()),
                         col.column("Transcript", DataSourceCreator.TRANSCRIPT_FIELD, type.stringType()).setWidth(250),
                         col.column("CDS", DataSourceCreator.HGVS_CODING_FIELD, type.stringType()),
                         col.column("AA", DataSourceCreator.HGVS_PROTEIN_FIELD, type.stringType()),
                         col.column("Effect", DataSourceCreator.EFFECT_FIELD, type.stringType()),
                         col.column("Cosmic", DataSourceCreator.COSMIC_FIELD, type.stringType()),
                         col.column("Read Count", DataSourceCreator.READ_COUNT_FIELD, type.stringType()))
                .setColumnStyle(stl.style()
                        .setFontSize(8)
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))
                .setColumnTitleStyle(stl.style()
                        .bold()
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                        .setBorder(stl.pen1Point())
                        .setBackgroundColor(Color.LIGHT_GRAY))
                .setDataSource(variantDataSource)
                .show()
                .print();
        // @formatter:on
    }

    private static class Item1Expression extends AbstractSimpleExpression<String> {
        @Override
        public String evaluate(ReportParameters reportParameters) {
            return reportParameters.getValue("item1");
        }
    }

    private static class Item2Expression extends AbstractSimpleExpression<String> {
        @Override
        public String evaluate(ReportParameters reportParameters) {
            return reportParameters.getValue("item2");
        }
    }
}
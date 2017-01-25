package com.hartwig.hmftools.patientreporter.output;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.col;
import static net.sf.dynamicreports.report.builder.DynamicReports.field;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberReport;
import com.hartwig.hmftools.patientreporter.variants.VariantReport;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

import net.sf.dynamicreports.report.base.expression.AbstractSimpleExpression;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.datasource.DRDataSource;
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
                .fields(
                        field("item1", String.class),
                        field("item2", String.class))
                .columns(col.componentColumn(listBuilder))
//                .columns(col.column("Gene", "gene", type.stringType()),
//                         col.column("Chr:Pos", "position", type.stringType()),
//                         col.column("Ref", "ref", type.stringType()),
//                         col.column("Alt", "alt", type.stringType()),
//                         col.column("HGVS Transcript", "transcript", type.stringType()).setWidth(250),
//                         col.column("HGVS Protein", "hgvsProtein", type.stringType()),
//                         col.column("Effect", "consequence", type.stringType()),
//                         col.column("Cosmic ID", "cosmicID", type.stringType()),
//                         col.column("Allele #", "alleleReadCount", type.integerType()),
//                         col.column("Total #", "totalReadCount", type.integerType()))
                .setColumnStyle(stl.style()
                        .setFontSize(8)
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER))
                .setColumnTitleStyle(stl.style()
                        .bold()
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                        .setBorder(stl.pen1Point())
                        .setBackgroundColor(Color.LIGHT_GRAY))
                .setDataSource(createDataSource())
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

    @NotNull
    private static JRDataSource createDataSource() {
        final DRDataSource dataSource = new DRDataSource("item1", "item2");
        dataSource.add("item1_1", "item1_2");
        dataSource.add("item2_1", "item2_2");

        return dataSource;
    }
}
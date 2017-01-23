package com.hartwig.hmftools.patientreporter.output;

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

import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;
import net.sf.dynamicreports.report.exception.DRException;

public class PDFWriterTest {

    private static final String HMF_LOGO = Resources.getResource("pdf/hartwig_logo.jpg").getPath();

    @Test
    @Ignore
    public void canGeneratePDFReport() throws DRException {
        final String sample = "CPCT11111111T";
        final List<VariantReport> variants = Lists.newArrayList(
                new VariantReport.Builder().gene("PIK3CA").position("3:178936082").ref("G").alt("A").transcript(
                        "ENST00000263967.3").hgvsCoding("c.1624G>A").hgvsProtein("p.Glu542Lys").consequence(
                        "missense_variant").cosmicID("COSM125369").alleleFrequency("77%").readDepth("151").build());

        final List<CopyNumberReport> copyNumbers = Lists.newArrayList();
        //        final int mutationalLoad = 0;
        //        final PatientReport patientReport = new PatientReport(sample, variants, copyNumbers, mutationalLoad);

        // @formatter:off
        report()
                .title(
                    cmp.horizontalList().add(
                        cmp.image(HMF_LOGO),
                        cmp.text("HMF Sequencing Report - " + sample).setStyle(stl.style()
                                .bold()
                                .setFontSize(12)
                                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE))))
                .columns(col.column("Gen", "gene", type.stringType()),
                         col.column("Positie", "position", type.stringType()),
                         col.column("Ref", "ref", type.stringType()),
                         col.column("Alt", "alt", type.stringType()))
                .setColumnTitleStyle(stl.style()
                        .bold()
                        .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                        .setBorder(stl.pen1Point())
                        .setBackgroundColor(Color.LIGHT_GRAY))
                .setDataSource(variants)
                .show()
                .print();
        // @formatter:on
    }
}
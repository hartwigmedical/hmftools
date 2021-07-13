package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.ImageUtil;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.jetbrains.annotations.NotNull;

public class QualityControlChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;

    public QualityControlChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Quality Control";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Quality Control").addStyle(ReportResources.chapterTitleStyle()));

        addKeyQC(document);
        addMetricsFlagstats(document);
        addPurpleQCPlots(document);
        addSageBQRPlots(document);
    }

    private void addKeyQC(@NotNull Document document) {
        document.add(new Paragraph("TODO: Add Key QC & Contamination").addStyle(ReportResources.tableContentStyle()));
    }

    private void addMetricsFlagstats(@NotNull Document document) {
        document.add(new Paragraph("TODO: Add Metrics / Flagstats").addStyle(ReportResources.tableContentStyle()));
    }

    private void addSageBQRPlots(@NotNull Document document) {
        document.add(new Paragraph("SAGE Reference Sample BQR plot").addStyle(ReportResources.tableTitleStyle()));
        Image refImage = ImageUtil.build(report.plots().sageReferenceBQRPlot());
        refImage.setMaxWidth(ReportResources.CONTENT_WIDTH);
        refImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(refImage);

        document.add(new Paragraph("SAGE Tumor Sample BQR plot").addStyle(ReportResources.tableTitleStyle()));
        Image tumorImage = ImageUtil.build(report.plots().sageTumorBQRPlot());
        tumorImage.setMaxWidth(ReportResources.CONTENT_WIDTH);
        tumorImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(tumorImage);
    }

    private void addPurpleQCPlots(@NotNull Document document) {
        document.add(new Paragraph("Purple QC plots").addStyle(ReportResources.tableTitleStyle()));

        long halfContentWidth = Math.round(ReportResources.CONTENT_WIDTH / 2D) - 2;
        Table table = new Table(2);
        table.addCell(TableUtil.createImageCell(ImageUtil.build(report.plots().purpleInputPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(TableUtil.createImageCell(ImageUtil.build(report.plots().purplePurityRangePlot()).setMaxWidth(halfContentWidth)));
        table.addCell(TableUtil.createImageCell(ImageUtil.build(report.plots().purpleCopyNumberPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(TableUtil.createImageCell(ImageUtil.build(report.plots().purpleVariantCopyNumberPlot())
                .setMaxWidth(halfContentWidth)));
        document.add(table);
    }
}

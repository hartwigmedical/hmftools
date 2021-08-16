package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.ImageUtil;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.jetbrains.annotations.NotNull;

public class QualityControlChapter implements ReportChapter {

    private static final DecimalFormat PERCENTAGE_FORMAT = ReportResources.decimalFormat("#'%'");

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

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));

        addKeyQC(document);
        addMetricsFlagstats(document);
        addPurpleQCPlots(document);
        addSageBQRPlots(document);
    }

    private void addKeyQC(@NotNull Document document) {
        Table table = TableUtil.createReportContentTable(new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("QC"), TableUtil.createHeaderCell("Amber Mean Depth"),
                        TableUtil.createHeaderCell("Contamination"), TableUtil.createHeaderCell("Uns. CN segments"),
                        TableUtil.createHeaderCell("Deleted Genes") });

        table.addCell(TableUtil.createContentCell(purpleQCString()));
        table.addCell(TableUtil.createContentCell(String.valueOf(report.purple().qc().amberMeanDepth())));
        table.addCell(TableUtil.createContentCell(PERCENTAGE_FORMAT.format(report.purple().qc().contamination() * 100)));
        table.addCell(TableUtil.createContentCell(String.valueOf(report.purple().qc().unsupportedCopyNumberSegments())));
        table.addCell(TableUtil.createContentCell(String.valueOf(report.purple().qc().deletedGenes())));

        document.add(TableUtil.createWrappingReportTable(table));
    }

    @NotNull
    private String purpleQCString() {
        StringJoiner joiner = new StringJoiner(", ");
        for (PurpleQCStatus status : report.purple().qc().status()) {
            joiner.add(status.toString());
        }
        return joiner.toString();
    }

    private void addMetricsFlagstats(@NotNull Document document) {
        document.add(new Paragraph("TODO: Add Metrics / Flagstats").addStyle(ReportResources.tableTitleStyle()));
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

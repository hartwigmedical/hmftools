package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.OrangeReport;
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

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class QualityControlChapter implements ReportChapter {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");
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
        addPurplePurityFitPlot(document);
        addFlagstats(document);
        addCoverageStats(document);
        addExcludedPercentages(document);
        addPurpleQCPlots(document);
        addSageBQRPlots(document);
    }

    private void addKeyQC(@NotNull Document document) {
        Table table = TableUtil.createReportContentTable(contentWidth(),
                new float[] { 2, 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("QC"), TableUtil.createHeaderCell("Fit Method"),
                        TableUtil.createHeaderCell("Mean Depth"), TableUtil.createHeaderCell("Contamination"),
                        TableUtil.createHeaderCell("Uns. CN segments"), TableUtil.createHeaderCell("Deleted Genes") });

        table.addCell(TableUtil.createContentCell(purpleQCString()));
        table.addCell(TableUtil.createContentCell(report.purple().fittedPurityMethod().toString()));
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

    private void addPurplePurityFitPlot(@NotNull Document document) {
        Image image = ImageUtil.build(report.plots().purplePurityRangePlot());
        image.setMaxWidth(contentWidth());
        image.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(image);
    }

    private void addFlagstats(@NotNull Document document) {
        Flagstat refFlagstat = report.refSample().flagstat();
        Flagstat tumorFlagstat = report.tumorSample().flagstat();

        Table flagstat = TableUtil.createReportContentTable(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell(Strings.EMPTY), TableUtil.createHeaderCell("Unique RC"),
                        TableUtil.createHeaderCell("Secondary RC"), TableUtil.createHeaderCell("Supplementary RC"),
                        TableUtil.createHeaderCell("Mapped Proportion") });

        flagstat.addCell(TableUtil.createContentCell("Ref Sample"));
        flagstat.addCell(TableUtil.createContentCell(String.valueOf(refFlagstat.uniqueReadCount())));
        flagstat.addCell(TableUtil.createContentCell(String.valueOf(refFlagstat.secondaryCount())));
        flagstat.addCell(TableUtil.createContentCell(String.valueOf(refFlagstat.supplementaryCount())));
        flagstat.addCell(TableUtil.createContentCell(PERCENTAGE_FORMAT.format(refFlagstat.mappedProportion() * 100)));

        flagstat.addCell(TableUtil.createContentCell("Tumor Sample"));
        flagstat.addCell(TableUtil.createContentCell(String.valueOf(tumorFlagstat.uniqueReadCount())));
        flagstat.addCell(TableUtil.createContentCell(String.valueOf(tumorFlagstat.secondaryCount())));
        flagstat.addCell(TableUtil.createContentCell(String.valueOf(tumorFlagstat.supplementaryCount())));
        flagstat.addCell(TableUtil.createContentCell(PERCENTAGE_FORMAT.format(tumorFlagstat.mappedProportion() * 100)));

        document.add(TableUtil.createWrappingReportTable(flagstat, "Flagstats"));
    }

    private void addCoverageStats(@NotNull Document document) {
        WGSMetrics refMetrics = report.refSample().metrics();
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Table coverage = TableUtil.createReportContentTable(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell(Strings.EMPTY), TableUtil.createHeaderCell("Mean Coverage"),
                        TableUtil.createHeaderCell("SD Coverage"), TableUtil.createHeaderCell("Median Coverage"),
                        TableUtil.createHeaderCell("Mad Coverage") });

        coverage.addCell(TableUtil.createContentCell("Ref Sample"));
        coverage.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(refMetrics.meanCoverage())));
        coverage.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(refMetrics.sdCoverage())));
        coverage.addCell(TableUtil.createContentCell(String.valueOf(refMetrics.medianCoverage())));
        coverage.addCell(TableUtil.createContentCell(String.valueOf(refMetrics.madCoverage())));

        coverage.addCell(TableUtil.createContentCell("Tumor Sample"));
        coverage.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(tumorMetrics.meanCoverage())));
        coverage.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(tumorMetrics.sdCoverage())));
        coverage.addCell(TableUtil.createContentCell(String.valueOf(tumorMetrics.medianCoverage())));
        coverage.addCell(TableUtil.createContentCell(String.valueOf(tumorMetrics.madCoverage())));

        document.add(TableUtil.createWrappingReportTable(coverage, "Coverage Stats"));
    }

    private void addExcludedPercentages(@NotNull Document document) {
        WGSMetrics refMetrics = report.refSample().metrics();
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Table percentages = TableUtil.createReportContentTable(contentWidth(),
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell(Strings.EMPTY), TableUtil.createHeaderCell("Adapter"),
                        TableUtil.createHeaderCell("BaseQ"), TableUtil.createHeaderCell("Capped"), TableUtil.createHeaderCell("Dupe"),
                        TableUtil.createHeaderCell("MapQ"), TableUtil.createHeaderCell("Overlap"), TableUtil.createHeaderCell("Unpaired"),
                        TableUtil.createHeaderCell("Total") });

        percentages.addCell(TableUtil.createContentCell("Ref Sample"));
        percentages.addCell(TableUtil.createContentCell(percent(refMetrics.pctExcAdapter())));
        percentages.addCell(TableUtil.createContentCell(percent(refMetrics.pctExcBaseQ())));
        percentages.addCell(TableUtil.createContentCell(percent(refMetrics.pctExcCapped())));
        percentages.addCell(TableUtil.createContentCell(percent(refMetrics.pctExcDupe())));
        percentages.addCell(TableUtil.createContentCell(percent(refMetrics.pctExcMapQ())));
        percentages.addCell(TableUtil.createContentCell(percent(refMetrics.pctExcOverlap())));
        percentages.addCell(TableUtil.createContentCell(percent(refMetrics.pctExcUnpaired())));
        percentages.addCell(TableUtil.createContentCell(percent(refMetrics.pctExcTotal())));

        percentages.addCell(TableUtil.createContentCell("Tumor Sample"));
        percentages.addCell(TableUtil.createContentCell(percent(tumorMetrics.pctExcAdapter())));
        percentages.addCell(TableUtil.createContentCell(percent(tumorMetrics.pctExcBaseQ())));
        percentages.addCell(TableUtil.createContentCell(percent(tumorMetrics.pctExcCapped())));
        percentages.addCell(TableUtil.createContentCell(percent(tumorMetrics.pctExcDupe())));
        percentages.addCell(TableUtil.createContentCell(percent(tumorMetrics.pctExcMapQ())));
        percentages.addCell(TableUtil.createContentCell(percent(tumorMetrics.pctExcOverlap())));
        percentages.addCell(TableUtil.createContentCell(percent(tumorMetrics.pctExcUnpaired())));
        percentages.addCell(TableUtil.createContentCell(percent(tumorMetrics.pctExcTotal())));

        document.add(TableUtil.createWrappingReportTable(percentages, "Excluded Percentages"));
    }

    @NotNull
    private static String percent(@Nullable Double value) {
        return value != null ? PERCENTAGE_FORMAT.format(value * 100) : ReportResources.NOT_AVAILABLE;
    }

    private void addPurpleQCPlots(@NotNull Document document) {
        document.add(new Paragraph("QC plots").addStyle(ReportResources.tableTitleStyle()));

        long halfContentWidth = Math.round(contentWidth() / 2D) - 2;
        Table table = new Table(2);
        table.addCell(TableUtil.createImageCell(ImageUtil.build(report.plots().purpleFinalCircosPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(TableUtil.createImageCell(ImageUtil.build(report.plots().purpleInputPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(TableUtil.createImageCell(ImageUtil.build(report.plots().purpleCopyNumberPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(TableUtil.createImageCell(ImageUtil.build(report.plots().purpleVariantCopyNumberPlot())
                .setMaxWidth(halfContentWidth)));
        document.add(table);
    }

    private void addSageBQRPlots(@NotNull Document document) {
        document.add(new Paragraph("Reference Sample BQR plot").addStyle(ReportResources.tableTitleStyle()));
        Image refImage = ImageUtil.build(report.plots().sageReferenceBQRPlot());
        refImage.setMaxWidth(contentWidth());
        refImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(refImage);

        document.add(new Paragraph("Tumor Sample BQR plot").addStyle(ReportResources.tableTitleStyle()));
        Image tumorImage = ImageUtil.build(report.plots().sageTumorBQRPlot());
        tumorImage.setMaxWidth(contentWidth());
        tumorImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(tumorImage);
    }
}

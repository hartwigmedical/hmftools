package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.CellUtil;
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
        Table table = TableUtil.createContent(contentWidth(),
                new float[] { 2, 1, 1, 1, 1, 1 },
                new Cell[] { CellUtil.createHeader("QC"), CellUtil.createHeader("Fit Method"), CellUtil.createHeader("Mean Depth"),
                        CellUtil.createHeader("Contamination"), CellUtil.createHeader("Uns. CN segments"),
                        CellUtil.createHeader("Deleted Genes") });

        table.addCell(CellUtil.createContent(purpleQCString()));
        table.addCell(CellUtil.createContent(report.purple().fittedPurityMethod().toString()));
        table.addCell(CellUtil.createContent(String.valueOf(report.purple().qc().amberMeanDepth())));
        table.addCell(CellUtil.createContent(PERCENTAGE_FORMAT.format(report.purple().qc().contamination() * 100)));
        table.addCell(CellUtil.createContent(String.valueOf(report.purple().qc().unsupportedCopyNumberSegments())));
        table.addCell(CellUtil.createContent(String.valueOf(report.purple().qc().deletedGenes())));

        document.add(TableUtil.createWrapping(table));
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

        Table flagstat = TableUtil.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { CellUtil.createHeader(Strings.EMPTY), CellUtil.createHeader("Unique RC"),
                        CellUtil.createHeader("Secondary RC"), CellUtil.createHeader("Supplementary RC"),
                        CellUtil.createHeader("Mapped Proportion") });

        flagstat.addCell(CellUtil.createContent("Ref Sample"));
        flagstat.addCell(CellUtil.createContent(String.valueOf(refFlagstat.uniqueReadCount())));
        flagstat.addCell(CellUtil.createContent(String.valueOf(refFlagstat.secondaryCount())));
        flagstat.addCell(CellUtil.createContent(String.valueOf(refFlagstat.supplementaryCount())));
        flagstat.addCell(CellUtil.createContent(PERCENTAGE_FORMAT.format(refFlagstat.mappedProportion() * 100)));

        flagstat.addCell(CellUtil.createContent("Tumor Sample"));
        flagstat.addCell(CellUtil.createContent(String.valueOf(tumorFlagstat.uniqueReadCount())));
        flagstat.addCell(CellUtil.createContent(String.valueOf(tumorFlagstat.secondaryCount())));
        flagstat.addCell(CellUtil.createContent(String.valueOf(tumorFlagstat.supplementaryCount())));
        flagstat.addCell(CellUtil.createContent(PERCENTAGE_FORMAT.format(tumorFlagstat.mappedProportion() * 100)));

        document.add(TableUtil.createWrapping(flagstat, "Flagstats"));
    }

    private void addCoverageStats(@NotNull Document document) {
        WGSMetrics refMetrics = report.refSample().metrics();
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Table coverage = TableUtil.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { CellUtil.createHeader(Strings.EMPTY), CellUtil.createHeader("Mean Coverage"),
                        CellUtil.createHeader("SD Coverage"), CellUtil.createHeader("Median Coverage"),
                        CellUtil.createHeader("Mad Coverage") });

        coverage.addCell(CellUtil.createContent("Ref Sample"));
        coverage.addCell(CellUtil.createContent(SINGLE_DIGIT.format(refMetrics.meanCoverage())));
        coverage.addCell(CellUtil.createContent(SINGLE_DIGIT.format(refMetrics.sdCoverage())));
        coverage.addCell(CellUtil.createContent(String.valueOf(refMetrics.medianCoverage())));
        coverage.addCell(CellUtil.createContent(String.valueOf(refMetrics.madCoverage())));

        coverage.addCell(CellUtil.createContent("Tumor Sample"));
        coverage.addCell(CellUtil.createContent(SINGLE_DIGIT.format(tumorMetrics.meanCoverage())));
        coverage.addCell(CellUtil.createContent(SINGLE_DIGIT.format(tumorMetrics.sdCoverage())));
        coverage.addCell(CellUtil.createContent(String.valueOf(tumorMetrics.medianCoverage())));
        coverage.addCell(CellUtil.createContent(String.valueOf(tumorMetrics.madCoverage())));

        document.add(TableUtil.createWrapping(coverage, "Coverage Stats"));
    }

    private void addExcludedPercentages(@NotNull Document document) {
        WGSMetrics refMetrics = report.refSample().metrics();
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Table percentages = TableUtil.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { CellUtil.createHeader(Strings.EMPTY), CellUtil.createHeader("Adapter"), CellUtil.createHeader("BaseQ"),
                        CellUtil.createHeader("Capped"), CellUtil.createHeader("Dupe"), CellUtil.createHeader("MapQ"),
                        CellUtil.createHeader("Overlap"), CellUtil.createHeader("Unpaired"), CellUtil.createHeader("Total") });

        percentages.addCell(CellUtil.createContent("Ref Sample"));
        percentages.addCell(CellUtil.createContent(percent(refMetrics.pctExcAdapter())));
        percentages.addCell(CellUtil.createContent(percent(refMetrics.pctExcBaseQ())));
        percentages.addCell(CellUtil.createContent(percent(refMetrics.pctExcCapped())));
        percentages.addCell(CellUtil.createContent(percent(refMetrics.pctExcDupe())));
        percentages.addCell(CellUtil.createContent(percent(refMetrics.pctExcMapQ())));
        percentages.addCell(CellUtil.createContent(percent(refMetrics.pctExcOverlap())));
        percentages.addCell(CellUtil.createContent(percent(refMetrics.pctExcUnpaired())));
        percentages.addCell(CellUtil.createContent(percent(refMetrics.pctExcTotal())));

        percentages.addCell(CellUtil.createContent("Tumor Sample"));
        percentages.addCell(CellUtil.createContent(percent(tumorMetrics.pctExcAdapter())));
        percentages.addCell(CellUtil.createContent(percent(tumorMetrics.pctExcBaseQ())));
        percentages.addCell(CellUtil.createContent(percent(tumorMetrics.pctExcCapped())));
        percentages.addCell(CellUtil.createContent(percent(tumorMetrics.pctExcDupe())));
        percentages.addCell(CellUtil.createContent(percent(tumorMetrics.pctExcMapQ())));
        percentages.addCell(CellUtil.createContent(percent(tumorMetrics.pctExcOverlap())));
        percentages.addCell(CellUtil.createContent(percent(tumorMetrics.pctExcUnpaired())));
        percentages.addCell(CellUtil.createContent(percent(tumorMetrics.pctExcTotal())));

        document.add(TableUtil.createWrapping(percentages, "Excluded Percentages"));
    }

    @NotNull
    private static String percent(@Nullable Double value) {
        return value != null ? PERCENTAGE_FORMAT.format(value * 100) : ReportResources.NOT_AVAILABLE;
    }

    private void addPurpleQCPlots(@NotNull Document document) {
        document.add(new Paragraph("QC plots").addStyle(ReportResources.tableTitleStyle()));

        long halfContentWidth = Math.round(contentWidth() / 2D) - 2;
        Table table = new Table(2);
        table.addCell(CellUtil.createImage(ImageUtil.build(report.plots().purpleFinalCircosPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(CellUtil.createImage(ImageUtil.build(report.plots().purpleInputPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(CellUtil.createImage(ImageUtil.build(report.plots().purpleCopyNumberPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(CellUtil.createImage(ImageUtil.build(report.plots().purpleVariantCopyNumberPlot()).setMaxWidth(halfContentWidth)));
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

package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Images;
import com.hartwig.hmftools.orange.report.util.Tables;
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
        Table table = Tables.createContent(contentWidth(),
                new float[] { 2, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("QC"), Cells.createHeader("Fit Method"), Cells.createHeader("Mean Depth"),
                        Cells.createHeader("Contamination"), Cells.createHeader("Uns. CN segments"),
                        Cells.createHeader("Deleted Genes") });

        table.addCell(Cells.createContent(purpleQCString()));
        table.addCell(Cells.createContent(report.purple().fittedPurityMethod().toString()));
        table.addCell(Cells.createContent(String.valueOf(report.purple().qc().amberMeanDepth())));
        table.addCell(Cells.createContent(PERCENTAGE_FORMAT.format(report.purple().qc().contamination() * 100)));
        table.addCell(Cells.createContent(String.valueOf(report.purple().qc().unsupportedCopyNumberSegments())));
        table.addCell(Cells.createContent(String.valueOf(report.purple().qc().deletedGenes())));

        document.add(Tables.createWrapping(table));
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
        Image image = Images.build(report.plots().purplePurityRangePlot());
        image.setMaxWidth(contentWidth());
        image.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(image);
    }

    private void addFlagstats(@NotNull Document document) {
        Flagstat refFlagstat = report.refSample().flagstat();
        Flagstat tumorFlagstat = report.tumorSample().flagstat();

        Table flagstat = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader(Strings.EMPTY), Cells.createHeader("Unique RC"),
                        Cells.createHeader("Secondary RC"), Cells.createHeader("Supplementary RC"),
                        Cells.createHeader("Mapped Proportion") });

        flagstat.addCell(Cells.createContent("Ref Sample"));
        flagstat.addCell(Cells.createContent(String.valueOf(refFlagstat.uniqueReadCount())));
        flagstat.addCell(Cells.createContent(String.valueOf(refFlagstat.secondaryCount())));
        flagstat.addCell(Cells.createContent(String.valueOf(refFlagstat.supplementaryCount())));
        flagstat.addCell(Cells.createContent(PERCENTAGE_FORMAT.format(refFlagstat.mappedProportion() * 100)));

        flagstat.addCell(Cells.createContent("Tumor Sample"));
        flagstat.addCell(Cells.createContent(String.valueOf(tumorFlagstat.uniqueReadCount())));
        flagstat.addCell(Cells.createContent(String.valueOf(tumorFlagstat.secondaryCount())));
        flagstat.addCell(Cells.createContent(String.valueOf(tumorFlagstat.supplementaryCount())));
        flagstat.addCell(Cells.createContent(PERCENTAGE_FORMAT.format(tumorFlagstat.mappedProportion() * 100)));

        document.add(Tables.createWrapping(flagstat, "Flagstats"));
    }

    private void addCoverageStats(@NotNull Document document) {
        WGSMetrics refMetrics = report.refSample().metrics();
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Table coverage = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader(Strings.EMPTY), Cells.createHeader("Mean Coverage"),
                        Cells.createHeader("SD Coverage"), Cells.createHeader("Median Coverage"),
                        Cells.createHeader("Mad Coverage") });

        coverage.addCell(Cells.createContent("Ref Sample"));
        coverage.addCell(Cells.createContent(SINGLE_DIGIT.format(refMetrics.meanCoverage())));
        coverage.addCell(Cells.createContent(SINGLE_DIGIT.format(refMetrics.sdCoverage())));
        coverage.addCell(Cells.createContent(String.valueOf(refMetrics.medianCoverage())));
        coverage.addCell(Cells.createContent(String.valueOf(refMetrics.madCoverage())));

        coverage.addCell(Cells.createContent("Tumor Sample"));
        coverage.addCell(Cells.createContent(SINGLE_DIGIT.format(tumorMetrics.meanCoverage())));
        coverage.addCell(Cells.createContent(SINGLE_DIGIT.format(tumorMetrics.sdCoverage())));
        coverage.addCell(Cells.createContent(String.valueOf(tumorMetrics.medianCoverage())));
        coverage.addCell(Cells.createContent(String.valueOf(tumorMetrics.madCoverage())));

        document.add(Tables.createWrapping(coverage, "Coverage Stats"));
    }

    private void addExcludedPercentages(@NotNull Document document) {
        WGSMetrics refMetrics = report.refSample().metrics();
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Table percentages = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader(Strings.EMPTY), Cells.createHeader("Adapter"), Cells.createHeader("BaseQ"),
                        Cells.createHeader("Capped"), Cells.createHeader("Dupe"), Cells.createHeader("MapQ"),
                        Cells.createHeader("Overlap"), Cells.createHeader("Unpaired"), Cells.createHeader("Total") });

        percentages.addCell(Cells.createContent("Ref Sample"));
        percentages.addCell(Cells.createContent(percent(refMetrics.pctExcAdapter())));
        percentages.addCell(Cells.createContent(percent(refMetrics.pctExcBaseQ())));
        percentages.addCell(Cells.createContent(percent(refMetrics.pctExcCapped())));
        percentages.addCell(Cells.createContent(percent(refMetrics.pctExcDupe())));
        percentages.addCell(Cells.createContent(percent(refMetrics.pctExcMapQ())));
        percentages.addCell(Cells.createContent(percent(refMetrics.pctExcOverlap())));
        percentages.addCell(Cells.createContent(percent(refMetrics.pctExcUnpaired())));
        percentages.addCell(Cells.createContent(percent(refMetrics.pctExcTotal())));

        percentages.addCell(Cells.createContent("Tumor Sample"));
        percentages.addCell(Cells.createContent(percent(tumorMetrics.pctExcAdapter())));
        percentages.addCell(Cells.createContent(percent(tumorMetrics.pctExcBaseQ())));
        percentages.addCell(Cells.createContent(percent(tumorMetrics.pctExcCapped())));
        percentages.addCell(Cells.createContent(percent(tumorMetrics.pctExcDupe())));
        percentages.addCell(Cells.createContent(percent(tumorMetrics.pctExcMapQ())));
        percentages.addCell(Cells.createContent(percent(tumorMetrics.pctExcOverlap())));
        percentages.addCell(Cells.createContent(percent(tumorMetrics.pctExcUnpaired())));
        percentages.addCell(Cells.createContent(percent(tumorMetrics.pctExcTotal())));

        document.add(Tables.createWrapping(percentages, "Excluded Percentages"));
    }

    @NotNull
    private static String percent(@Nullable Double value) {
        return value != null ? PERCENTAGE_FORMAT.format(value * 100) : ReportResources.NOT_AVAILABLE;
    }

    private void addPurpleQCPlots(@NotNull Document document) {
        document.add(new Paragraph("QC plots").addStyle(ReportResources.tableTitleStyle()));

        long halfContentWidth = Math.round(contentWidth() / 2D) - 2;
        Table table = new Table(2);
        table.addCell(Cells.createImage(Images.build(report.plots().purpleFinalCircosPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(Cells.createImage(Images.build(report.plots().purpleInputPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(Cells.createImage(Images.build(report.plots().purpleCopyNumberPlot()).setMaxWidth(halfContentWidth)));
        table.addCell(Cells.createImage(Images.build(report.plots().purpleVariantCopyNumberPlot()).setMaxWidth(halfContentWidth)));
        document.add(table);
    }

    private void addSageBQRPlots(@NotNull Document document) {
        document.add(new Paragraph("Reference Sample BQR plot").addStyle(ReportResources.tableTitleStyle()));
        Image refImage = Images.build(report.plots().sageReferenceBQRPlot());
        refImage.setMaxWidth(contentWidth());
        refImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(refImage);

        document.add(new Paragraph("Tumor Sample BQR plot").addStyle(ReportResources.tableTitleStyle()));
        Image tumorImage = Images.build(report.plots().sageTumorBQRPlot());
        tumorImage.setMaxWidth(contentWidth());
        tumorImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(tumorImage);
    }
}

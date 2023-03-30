package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentage;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.StringJoiner;

import com.hartwig.hmftools.datamodel.flagstat.Flagstat;
import com.hartwig.hmftools.datamodel.metrics.WGSMetrics;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
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

    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final PlotPathResolver plotPathResolver;

    public QualityControlChapter(@NotNull final OrangeRecord report, @NotNull final PlotPathResolver plotPathResolver) {
        this.report = report;
        this.plotPathResolver = plotPathResolver;
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
                new float[] { 2, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("QC"), Cells.createHeader("Ref Genome"), Cells.createHeader("Fit Method"),
                        Cells.createHeader("Mean Depth"), Cells.createHeader("Contamination"), Cells.createHeader("Uns. Segments"),
                        Cells.createHeader("Deleted Genes") });

        table.addCell(Cells.createContent(purpleQCString()));
        table.addCell(Cells.createContent(report.refGenomeVersion().toString()));
        table.addCell(Cells.createContent(report.purple().fit().fittedPurityMethod().toString()));
        table.addCell(Cells.createContent(String.valueOf(report.purple().fit().qc().amberMeanDepth())));
        table.addCell(Cells.createContent(formatPercentage(report.purple().fit().qc().contamination())));
        table.addCell(Cells.createContent(String.valueOf(report.purple().fit().qc().unsupportedCopyNumberSegments())));
        table.addCell(Cells.createContent(String.valueOf(report.purple().fit().qc().deletedGenes())));

        document.add(Tables.createWrapping(table));
    }

    @NotNull
    private String purpleQCString() {
        StringJoiner joiner = new StringJoiner(", ");
        for (PurpleQCStatus status : report.purple().fit().qc().status()) {
            joiner.add(status.toString());
        }
        return joiner.toString();
    }

    private void addPurplePurityFitPlot(@NotNull Document document) {
        Image image = Images.build(plotPathResolver.resolve(report.plots().purplePurityRangePlot()));
        image.setMaxWidth(contentWidth());
        image.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(image);
    }

    private void addFlagstats(@NotNull Document document) {
        Flagstat refFlagstat = report.refSample() != null ? report.refSample().flagstat() : null;
        Flagstat tumorFlagstat = report.tumorSample().flagstat();

        Table flagstat = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader(Strings.EMPTY), Cells.createHeader("Unique RC"), Cells.createHeader("Secondary RC"),
                        Cells.createHeader("Supplementary RC"), Cells.createHeader("Mapped Proportion") });

        if (refFlagstat != null) {
            flagstat.addCell(Cells.createContent("Ref Sample"));
            flagstat.addCell(Cells.createContent(String.valueOf(refFlagstat.uniqueReadCount())));
            flagstat.addCell(Cells.createContent(String.valueOf(refFlagstat.secondaryCount())));
            flagstat.addCell(Cells.createContent(String.valueOf(refFlagstat.supplementaryCount())));
            flagstat.addCell(Cells.createContent(formatPercentage(refFlagstat.mappedProportion())));
        }

        flagstat.addCell(Cells.createContent("Tumor Sample"));
        flagstat.addCell(Cells.createContent(String.valueOf(tumorFlagstat.uniqueReadCount())));
        flagstat.addCell(Cells.createContent(String.valueOf(tumorFlagstat.secondaryCount())));
        flagstat.addCell(Cells.createContent(String.valueOf(tumorFlagstat.supplementaryCount())));
        flagstat.addCell(Cells.createContent(formatPercentage(tumorFlagstat.mappedProportion())));

        document.add(Tables.createWrapping(flagstat, "Flagstats"));
    }

    private void addCoverageStats(@NotNull Document document) {
        WGSMetrics refMetrics = report.refSample() != null ? report.refSample().metrics() : null;
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Table coverage = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader(Strings.EMPTY), Cells.createHeader("Mean Coverage"), Cells.createHeader("SD Coverage"),
                        Cells.createHeader("Median Coverage"), Cells.createHeader("Mad Coverage") });

        if (refMetrics != null) {
            coverage.addCell(Cells.createContent("Ref Sample"));
            coverage.addCell(Cells.createContent(formatSingleDigitDecimal(refMetrics.meanCoverage())));
            coverage.addCell(Cells.createContent(formatSingleDigitDecimal(refMetrics.sdCoverage())));
            coverage.addCell(Cells.createContent(String.valueOf(refMetrics.medianCoverage())));
            coverage.addCell(Cells.createContent(String.valueOf(refMetrics.madCoverage())));
        }

        coverage.addCell(Cells.createContent("Tumor Sample"));
        coverage.addCell(Cells.createContent(formatSingleDigitDecimal(tumorMetrics.meanCoverage())));
        coverage.addCell(Cells.createContent(formatSingleDigitDecimal(tumorMetrics.sdCoverage())));
        coverage.addCell(Cells.createContent(String.valueOf(tumorMetrics.medianCoverage())));
        coverage.addCell(Cells.createContent(String.valueOf(tumorMetrics.madCoverage())));

        document.add(Tables.createWrapping(coverage, "Coverage Stats"));
    }

    private void addExcludedPercentages(@NotNull Document document) {
        WGSMetrics refMetrics = report.refSample() != null ? report.refSample().metrics() : null;
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Table percentages = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader(Strings.EMPTY), Cells.createHeader("Adapter"), Cells.createHeader("BaseQ"),
                        Cells.createHeader("Capped"), Cells.createHeader("Dupe"), Cells.createHeader("MapQ"), Cells.createHeader("Overlap"),
                        Cells.createHeader("Unpaired"), Cells.createHeader("Total") });

        if (refMetrics != null) {
            percentages.addCell(Cells.createContent("Ref Sample"));
            percentages.addCell(Cells.createContent(percent(refMetrics.pctExcAdapter())));
            percentages.addCell(Cells.createContent(percent(refMetrics.pctExcBaseQ())));
            percentages.addCell(Cells.createContent(percent(refMetrics.pctExcCapped())));
            percentages.addCell(Cells.createContent(percent(refMetrics.pctExcDupe())));
            percentages.addCell(Cells.createContent(percent(refMetrics.pctExcMapQ())));
            percentages.addCell(Cells.createContent(percent(refMetrics.pctExcOverlap())));
            percentages.addCell(Cells.createContent(percent(refMetrics.pctExcUnpaired())));
            percentages.addCell(Cells.createContent(percent(refMetrics.pctExcTotal())));
        }

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
        return value != null ? formatPercentage(value) : ReportResources.NOT_AVAILABLE;
    }

    private void addPurpleQCPlots(@NotNull Document document) {
        document.add(new Paragraph("QC plots").addStyle(ReportResources.tableTitleStyle()));

        long halfContentWidth = Math.round(contentWidth() / 2D) - 2;
        Table table = new Table(2);
        table.addCell(Cells.createImage(Images.build(plotPathResolver.resolve(report.plots().purpleFinalCircosPlot()))
                .setMaxWidth(halfContentWidth)));
        table.addCell(Cells.createImage(Images.build(plotPathResolver.resolve(report.plots().purpleInputPlot()))
                .setMaxWidth(halfContentWidth)));
        table.addCell(Cells.createImage(Images.build(plotPathResolver.resolve(report.plots().purpleCopyNumberPlot()))
                .setMaxWidth(halfContentWidth)));
        table.addCell(Cells.createImage(Images.build(plotPathResolver.resolve(report.plots().purpleVariantCopyNumberPlot()))
                .setMaxWidth(halfContentWidth)));
        document.add(table);
    }

    private void addSageBQRPlots(@NotNull Document document) {
        String sageReferenceBQRPlot = report.plots().sageReferenceBQRPlot();
        if (sageReferenceBQRPlot != null) {
            document.add(new Paragraph("Reference Sample BQR plot").addStyle(ReportResources.tableTitleStyle()));
            Image refImage = Images.build(plotPathResolver.resolve(sageReferenceBQRPlot));
            refImage.setMaxWidth(contentWidth());
            refImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
            document.add(refImage);
        }

        document.add(new Paragraph("Tumor Sample BQR plot").addStyle(ReportResources.tableTitleStyle()));
        Image tumorImage = Images.build(plotPathResolver.resolve(report.plots().sageTumorBQRPlot()));
        tumorImage.setMaxWidth(contentWidth());
        tumorImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(tumorImage);
    }
}

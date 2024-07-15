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

public class QualityControlChapter implements ReportChapter
{
    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final PlotPathResolver plotPathResolver;
    @NotNull
    private final ReportResources reportResources;

    public QualityControlChapter(@NotNull final OrangeRecord report, @NotNull final PlotPathResolver plotPathResolver,
            @NotNull final ReportResources reportResources)
    {
        this.report = report;
        this.plotPathResolver = plotPathResolver;
        this.reportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Quality Control";
    }

    @NotNull
    @Override
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document)
    {
        document.add(new Paragraph(name()).addStyle(reportResources.chapterTitleStyle()));

        addKeyQC(document);
        addPurplePurityFitPlot(document);
        addFlagstats(document);
        addCoverageStats(document);
        addExcludedPercentages(document);
        addPurpleQCPlots(document);
        addSageBQRPlots(document);
    }

    private void addKeyQC(@NotNull Document document)
    {
        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(contentWidth(),
                new float[] { 15, 10, 10, 10, 10, 12, 10 },
                new Cell[] { cells.createHeader("QC"), cells.createHeader("Ref Genome"), cells.createHeader("Fit Method"),
                        cells.createHeader("Mean Depth"), cells.createHeader("Contamination"), cells.createHeader("Uns. Segments (%)"),
                        cells.createHeader("Deleted Genes") });

        table.addCell(cells.createContent(purpleQCStatusString()));
        table.addCell(cells.createContent(report.refGenomeVersion().toString()));
        table.addCell(cells.createContent(report.purple().fit().fittedPurityMethod().toString()));
        table.addCell(cells.createContent(String.valueOf(report.purple().fit().qc().amberMeanDepth())));
        table.addCell(cells.createContent(formatPercentage(report.purple().fit().qc().contamination())));
        table.addCell(cells.createContent(purpleQCSegmentsString()));
        table.addCell(cells.createContent(String.valueOf(report.purple().fit().qc().deletedGenes())));

        document.add(new Tables(reportResources).createWrapping(table));
    }

    @NotNull
    private String purpleQCStatusString()
    {
        StringJoiner joiner = new StringJoiner(", ");
        for(PurpleQCStatus status : report.purple().fit().qc().status())
        {
            joiner.add(status.toString());
        }
        return joiner.toString();
    }

    @NotNull
    private String purpleQCSegmentsString()
    {
        int unsupportedSegments = report.purple().fit().qc().unsupportedCopyNumberSegments();
        return unsupportedSegments + " (" + percent(
                (double) unsupportedSegments / report.purple().fit().qc().totalCopyNumberSegments() * 100) + ")";
    }

    private void addPurplePurityFitPlot(@NotNull Document document)
    {
        Image image = Images.build(plotPathResolver.resolve(report.plots().purplePurityRangePlot()));
        image.setMaxWidth(contentWidth());
        image.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(image);
    }

    private void addFlagstats(@NotNull Document document)
    {
        Flagstat refFlagstat = report.refSample() != null ? report.refSample().flagstat() : null;
        Flagstat tumorFlagstat = report.tumorSample().flagstat();

        Cells cells = new Cells(reportResources);
        Table flagstat = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { cells.createHeader(Strings.EMPTY), cells.createHeader("Unique RC"), cells.createHeader("Secondary RC"),
                        cells.createHeader("Supplementary RC"), cells.createHeader("Mapped Proportion") });

        if(refFlagstat != null)
        {
            flagstat.addCell(cells.createContent("Ref Sample"));
            flagstat.addCell(cells.createContent(String.valueOf(refFlagstat.uniqueReadCount())));
            flagstat.addCell(cells.createContent(String.valueOf(refFlagstat.secondaryCount())));
            flagstat.addCell(cells.createContent(String.valueOf(refFlagstat.supplementaryCount())));
            flagstat.addCell(cells.createContent(formatPercentage(refFlagstat.mappedProportion())));
        }

        flagstat.addCell(cells.createContent("Tumor Sample"));
        flagstat.addCell(cells.createContent(String.valueOf(tumorFlagstat.uniqueReadCount())));
        flagstat.addCell(cells.createContent(String.valueOf(tumorFlagstat.secondaryCount())));
        flagstat.addCell(cells.createContent(String.valueOf(tumorFlagstat.supplementaryCount())));
        flagstat.addCell(cells.createContent(formatPercentage(tumorFlagstat.mappedProportion())));

        document.add(new Tables(reportResources).createWrapping(flagstat, "Flagstats"));
    }

    private void addCoverageStats(@NotNull Document document)
    {
        WGSMetrics refMetrics = report.refSample() != null ? report.refSample().metrics() : null;
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Cells cells = new Cells(reportResources);
        Table coverage = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { cells.createHeader(Strings.EMPTY), cells.createHeader("Mean Coverage"), cells.createHeader("SD Coverage"),
                        cells.createHeader("Median Coverage"), cells.createHeader("Mad Coverage") });

        if(refMetrics != null)
        {
            coverage.addCell(cells.createContent("Ref Sample"));
            coverage.addCell(cells.createContent(formatSingleDigitDecimal(refMetrics.meanCoverage())));
            coverage.addCell(cells.createContent(formatSingleDigitDecimal(refMetrics.sdCoverage())));
            coverage.addCell(cells.createContent(String.valueOf(refMetrics.medianCoverage())));
            coverage.addCell(cells.createContent(String.valueOf(refMetrics.madCoverage())));
        }

        coverage.addCell(cells.createContent("Tumor Sample"));
        coverage.addCell(cells.createContent(formatSingleDigitDecimal(tumorMetrics.meanCoverage())));
        coverage.addCell(cells.createContent(formatSingleDigitDecimal(tumorMetrics.sdCoverage())));
        coverage.addCell(cells.createContent(String.valueOf(tumorMetrics.medianCoverage())));
        coverage.addCell(cells.createContent(String.valueOf(tumorMetrics.madCoverage())));

        document.add(new Tables(reportResources).createWrapping(coverage, "Coverage Stats"));
    }

    private void addExcludedPercentages(@NotNull Document document)
    {
        WGSMetrics refMetrics = report.refSample() != null ? report.refSample().metrics() : null;
        WGSMetrics tumorMetrics = report.tumorSample().metrics();

        Cells cells = new Cells(reportResources);
        Table percentages = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { cells.createHeader(Strings.EMPTY), cells.createHeader("Adapter"), cells.createHeader("BaseQ"),
                        cells.createHeader("Capped"), cells.createHeader("Dupe"), cells.createHeader("MapQ"), cells.createHeader("Overlap"),
                        cells.createHeader("Unpaired"), cells.createHeader("Total") });

        if(refMetrics != null)
        {
            percentages.addCell(cells.createContent("Ref Sample"));
            percentages.addCell(cells.createContent(percent(refMetrics.pctExcAdapter())));
            percentages.addCell(cells.createContent(percent(refMetrics.pctExcBaseQ())));
            percentages.addCell(cells.createContent(percent(refMetrics.pctExcCapped())));
            percentages.addCell(cells.createContent(percent(refMetrics.pctExcDupe())));
            percentages.addCell(cells.createContent(percent(refMetrics.pctExcMapQ())));
            percentages.addCell(cells.createContent(percent(refMetrics.pctExcOverlap())));
            percentages.addCell(cells.createContent(percent(refMetrics.pctExcUnpaired())));
            percentages.addCell(cells.createContent(percent(refMetrics.pctExcTotal())));
        }

        percentages.addCell(cells.createContent("Tumor Sample"));
        percentages.addCell(cells.createContent(percent(tumorMetrics.pctExcAdapter())));
        percentages.addCell(cells.createContent(percent(tumorMetrics.pctExcBaseQ())));
        percentages.addCell(cells.createContent(percent(tumorMetrics.pctExcCapped())));
        percentages.addCell(cells.createContent(percent(tumorMetrics.pctExcDupe())));
        percentages.addCell(cells.createContent(percent(tumorMetrics.pctExcMapQ())));
        percentages.addCell(cells.createContent(percent(tumorMetrics.pctExcOverlap())));
        percentages.addCell(cells.createContent(percent(tumorMetrics.pctExcUnpaired())));
        percentages.addCell(cells.createContent(percent(tumorMetrics.pctExcTotal())));

        document.add(new Tables(reportResources).createWrapping(percentages, "Excluded Percentages"));
    }

    @NotNull
    private static String percent(@Nullable Double value)
    {
        return value != null ? formatPercentage(value) : ReportResources.NOT_AVAILABLE;
    }

    private void addPurpleQCPlots(@NotNull Document document)
    {
        document.add(new Paragraph("QC plots").addStyle(reportResources.tableTitleStyle()));

        long halfContentWidth = Math.round(contentWidth() / 2D) - 2;
        Table table = new Table(2);
        Cells cells = new Cells(reportResources);
        table.addCell(cells.createImage(Images.build(plotPathResolver.resolve(report.plots().purpleFinalCircosPlot()))
                .setMaxWidth(halfContentWidth)));
        table.addCell(cells.createImage(Images.build(plotPathResolver.resolve(report.plots().purpleInputPlot()))
                .setMaxWidth(halfContentWidth)));
        table.addCell(cells.createImage(Images.build(plotPathResolver.resolve(report.plots().purpleCopyNumberPlot()))
                .setMaxWidth(halfContentWidth)));
        table.addCell(cells.createImage(Images.build(plotPathResolver.resolve(report.plots().purpleVariantCopyNumberPlot()))
                .setMaxWidth(halfContentWidth)));
        document.add(table);
    }

    private void addSageBQRPlots(@NotNull Document document)
    {
        String sageReferenceBQRPlot = report.plots().sageReferenceBQRPlot();
        if(sageReferenceBQRPlot != null)
        {
            document.add(new Paragraph("Reference Sample BQR plot").addStyle(reportResources.tableTitleStyle()));
            Image refImage = Images.build(plotPathResolver.resolve(sageReferenceBQRPlot));
            refImage.setMaxWidth(contentWidth());
            refImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
            document.add(refImage);
        }

        document.add(new Paragraph("Tumor Sample BQR plot").addStyle(reportResources.tableTitleStyle()));
        Image tumorImage = Images.build(plotPathResolver.resolve(report.plots().sageTumorBQRPlot()));
        tumorImage.setMaxWidth(contentWidth());
        tumorImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(tumorImage);
    }
}

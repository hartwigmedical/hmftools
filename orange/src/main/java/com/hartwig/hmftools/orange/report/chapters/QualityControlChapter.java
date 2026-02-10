package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentage;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTwoDigitDecimal;

import java.util.StringJoiner;

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
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class QualityControlChapter implements ReportChapter
{
    private final OrangeRecord mReport;
    private final PlotPathResolver mPlotPathResolver;
    private final ReportResources mReportResources;

    public QualityControlChapter(final OrangeRecord report, final PlotPathResolver plotPathResolver,
            final ReportResources reportResources)
    {
        mReport = report;
        mPlotPathResolver = plotPathResolver;
        mReportResources = reportResources;
    }

    @Override
    public String name()
    {
        return "Quality Control";
    }

    @Override
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        addKeyQC(document);
        addPurplePurityFitPlot(document);

        if(!mReport.tumorOnlyMode())
        {
            addTumorStats(document);
        }
    }

    private void addKeyQC(final Document document)
    {
        Cells cells = new Cells(mReportResources);
        Table table = Tables.createContent(contentWidth(),
                new float[] { 15, 10, 10, 10, 10, 12, 10 },
                new Cell[] { cells.createHeader("QC"), cells.createHeader("Ref Genome"), cells.createHeader("Fit Method"),
                        cells.createHeader("Mean Depth"), cells.createHeader("Contamination"), cells.createHeader("Uns. Segments (%)"),
                        cells.createHeader("Deleted Genes") });

        table.addCell(cells.createContent(purpleQCStatusString()));
        table.addCell(cells.createContent(mReport.refGenomeVersion().toString()));
        table.addCell(cells.createContent(mReport.purple().fit().fittedPurityMethod().toString()));
        table.addCell(cells.createContent(String.valueOf(mReport.purple().fit().qc().amberMeanDepth())));
        table.addCell(cells.createContent(formatPercentage(mReport.purple().fit().qc().contamination())));
        table.addCell(cells.createContent(purpleQCSegmentsString()));
        table.addCell(cells.createContent(String.valueOf(mReport.purple().fit().qc().deletedGenes())));

        document.add(new Tables(mReportResources).createWrapping(table));
    }

    @NotNull
    private String purpleQCStatusString()
    {
        StringJoiner joiner = new StringJoiner(", ");
        for(PurpleQCStatus status : mReport.purple().fit().qc().status())
        {
            joiner.add(status.toString());
        }
        return joiner.toString();
    }

    @NotNull
    private String purpleQCSegmentsString()
    {
        int unsupportedSegments = mReport.purple().fit().qc().unsupportedCopyNumberSegments();
        return unsupportedSegments + " (" + percent(
                (double) unsupportedSegments / mReport.purple().fit().qc().totalCopyNumberSegments()) + ")";
    }

    private void addPurplePurityFitPlot(final Document document)
    {
        Image image = Images.build(mPlotPathResolver.resolve(mReport.plots().purplePurityRangePlot()));
        image.setMaxWidth(contentWidth());
        image.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(image);
    }

    @NotNull
    private static String percent(@Nullable Double value)
    {
        return value != null ? formatPercentage(value) : ReportResources.NOT_AVAILABLE;
    }


    private void addTumorStats(final Document document)
    {
        Cells cells = new Cells(mReportResources);
        Table tumorStats = new Table(UnitValue.createPercentArray(new float[] { 3, 1 })).setWidth(contentWidth());

        tumorStats.addCell(cells.createContent("Tumor maximum diploid proportion"));
        tumorStats.addCell(cells.createContent(formatTwoDigitDecimal(mReport.purple().tumorStats().maxDiploidProportion())));

        tumorStats.addCell(cells.createContent("Number of hotspot mutations"));
        tumorStats.addCell(cells.createContent(String.valueOf(mReport.purple().tumorStats().hotspotMutationCount())));

        tumorStats.addCell(cells.createContent("Sum of small variant allele read counts"));
        tumorStats.addCell(cells.createContent(String.valueOf(mReport.purple().tumorStats().smallVariantAlleleReadCount())));

        /* SV and segment information no longer available
        tumorStats.addCell(cells.createContent("Number of hotspot structural variants"));
        tumorStats.addCell(cells.createContent(String.valueOf(report.purple().tumorStats().hotspotStructuralVariantCount())));

        tumorStats.addCell(cells.createContent("Sum of structural variant tumor fragment counts (excluding single breakends)"));
        tumorStats.addCell(cells.createContent(String.valueOf(report.purple().tumorStats().structuralVariantTumorFragmentCount())));

        tumorStats.addCell(cells.createContent("Sum of B-allele frequency points in germline diploid regions with tumor ratio < 0.8 OR > 1.2"));
        tumorStats.addCell(cells.createContent(String.valueOf(report.purple().tumorStats().bafCount())));
        */

        document.add(new Tables(mReportResources).createWrapping(tumorStats, "Tumor Detection Statistics"));
    }
}

package com.hartwig.hmftools.orange.report.chapters;

import static java.lang.Math.round;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Images;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.apache.logging.log4j.util.Strings;

public class PurplePlotsChapter implements ReportChapter
{
    private final OrangeRecord mReport;
    private final PlotPathResolver mPlotPathResolver;
    private final ReportResources mReportResources;

    public PurplePlotsChapter(final OrangeRecord report, final PlotPathResolver plotPathResolver, final ReportResources reportResources)
    {
        mReport = report;
        mPlotPathResolver = plotPathResolver;
        mReportResources = reportResources;
    }

    @Override
    public String name()
    {
        return "Purity and Ploidy";
    }

    @Override
    public PageSize pageSize()
    {
        return PageSize.A4.rotate();
    }

    @Override
    public void render(final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        if(QcStatusInterpretation.hasPurpleFail(mReport.purple().fit().qc()))
        {
            mReportResources.addQcFailNotice(document);
            return;
        }

        addPurlePlots(document);
    }

    private static final int PLOT_IMAGE_HEIGHT = 225;

    private void addPurlePlots(final Document document)
    {
        Table table = new Table(3);
        Cells cells = new Cells(mReportResources);

        // layout:
        // row 1: input circos              CN PDF                  somatic clonality
        // row 2: ploidy/purity range       minor allele PDF        somatic rainfall

        float baseWidth = round((contentWidth() / 3D) - 2);
        float squareWidth = round(baseWidth * 0.8);
        float rectangeWidth = round(baseWidth * 1.2);

        addTableImage(table, cells, mReport.plots().purpleInputCircosPlot(), squareWidth);
        addTableImage(table, cells, mReport.plots().purpleCopyNumberPlot(), squareWidth);
        addTableImage(table, cells, mReport.plots().purpleClonalityPlot(), rectangeWidth);
        addTableImage(table, cells, mReport.plots().purplePurityRangePlot(), squareWidth);
        addTableImage(table, cells, mReport.plots().purpleMinorAlleleMapPlot(), squareWidth);
        addTableImage(table, cells, mReport.plots().purpleRainfallPlot(), rectangeWidth);

        // table.setBor
        // .setBorder(Border.NO_BORDER)

        document.add(table);
    }

    private void addTableImage(final Table table, final Cells cells, final String plotFilename, final float width)
    {
        Image image = Images.build(mPlotPathResolver.resolve(plotFilename));
        // float imageHeight = round((contentHeight() / 2D) - 2);
        // image.setMaxHeight(imageHeight);
        image.setMaxHeight(PLOT_IMAGE_HEIGHT);
        image.setHorizontalAlignment(HorizontalAlignment.CENTER);
        image.setMaxWidth(width);
        table.addCell(cells.createImage(image));
    }
}

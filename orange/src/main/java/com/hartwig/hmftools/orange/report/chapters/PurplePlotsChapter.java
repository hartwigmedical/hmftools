package com.hartwig.hmftools.orange.report.chapters;

import static java.lang.Math.round;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.PurplePlotsData;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Images;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.jetbrains.annotations.NotNull;

public class PurplePlotsChapter implements ReportChapter
{
    private final PurplePlotsData mData;
    private final ReportResources mReportResources;

    public PurplePlotsChapter(final PurplePlotsData data, final ReportResources reportResources)
    {
        mData = data;
        mReportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Purity and Ploidy";
    }

    @NotNull
    @Override
    public PageSize pageSize()
    {
        return PageSize.A4.rotate();
    }

    @Override
    public void render(@NotNull final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        if(mData.hasPurpleFail)
        {
            mReportResources.addQcFailNotice(document);
            return;
        }

        addPurplePlots(document);
    }

    private static final int PLOT_IMAGE_HEIGHT = 225;

    private void addPurplePlots(final Document document)
    {
        Table table = new Table(3);
        Cells cells = new Cells(mReportResources);

        // layout:
        // row 1: input circos              CN PDF                  somatic clonality
        // row 2: ploidy/purity range       minor allele PDF        somatic rainfall

        float baseWidth = round((contentWidth() / 3D) - 2);
        float squareWidth = round(baseWidth * 0.8);
        float rectangeWidth = round(baseWidth * 1.2);

        addTableImage(table, cells, mData.purpleInputCircosPlotPath, squareWidth);
        addTableImage(table, cells, mData.purpleCopyNumberPlotPath, squareWidth);
        addTableImage(table, cells, mData.purpleClonalityPlotPath, rectangeWidth);
        addTableImage(table, cells, mData.purplePurityRangePlotPath, squareWidth);
        addTableImage(table, cells, mData.purpleMinorAlleleMapPlotPath, squareWidth);
        addTableImage(table, cells, mData.purpleRainfallPlotPath, rectangeWidth);

        document.add(table);
    }

    private void addTableImage(final Table table, final Cells cells, final String plotPath, final float width)
    {
        Image image = Images.build(plotPath);
        image.setMaxHeight(PLOT_IMAGE_HEIGHT);
        image.setHorizontalAlignment(HorizontalAlignment.CENTER);
        image.setMaxWidth(width);
        table.addCell(cells.createImage(image));
    }
}

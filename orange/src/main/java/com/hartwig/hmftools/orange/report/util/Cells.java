package com.hartwig.hmftools.orange.report.util;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.borders.SolidBorder;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.IBlockElement;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.VerticalAlignment;

import org.jetbrains.annotations.NotNull;

public class Cells
{
    @NotNull
    private final ReportResources reportResources;

    public Cells(@NotNull ReportResources reportResources)
    {
        this.reportResources = reportResources;
    }

    @NotNull
    public Cell createHeader(@NotNull String text)
    {
        Cell cell = new Cell(1, 1);
        cell.setHeight(23); // Set fixed height to create consistent spacing between table title and header
        cell.setBorder(Border.NO_BORDER);
        cell.setVerticalAlignment(VerticalAlignment.BOTTOM);
        cell.addStyle(reportResources.tableHeaderStyle());
        cell.add(new Paragraph(text.toUpperCase()));
        return cell;
    }

    @NotNull
    public Cell createSpanningEntry(@NotNull Table table, @NotNull String text)
    {
        Cell cell = new Cell(1, table.getNumberOfColumns());
        cell.add(new Paragraph(text));
        cell.setBorder(Border.NO_BORDER);
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(reportResources.tableContentStyle());
        return cell;
    }

    @NotNull
    public Cell createSpanningWarning(@NotNull Table table, @NotNull String text)
    {
        Cell cell = new Cell(1, table.getNumberOfColumns());
        cell.add(new Paragraph(text));
        cell.setBorder(Border.NO_BORDER);
        cell.addStyle(reportResources.qcWarningStyle());
        return cell;
    }

    @NotNull
    public Cell createUrl(@NotNull String text, @NotNull String url)
    {
        Cell cell = createBorderlessBase();
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(reportResources.urlStyle());
        cell.add(new Paragraph(text));
        cell.setAction(PdfAction.createURI(url));
        return cell;
    }

    @NotNull
    public Cell createContent(@NotNull String text)
    {
        return createContent(new Paragraph(text));
    }

    @NotNull
    public Cell createContent(@NotNull IBlockElement element)
    {
        Cell cell = createBorderlessBase();
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(reportResources.tableContentStyle());
        cell.add(element);
        return cell;
    }

    @NotNull
    public Cell createImage(@NotNull Image image)
    {
        Cell cell = createBorderlessBase();
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(reportResources.tableContentStyle());
        cell.add(image);
        return cell;
    }

    @NotNull
    public Cell createKey(@NotNull String text)
    {
        Cell cell = createBorderlessBase();
        cell.addStyle(reportResources.keyStyle());
        cell.add(new Paragraph(text));
        return cell;
    }

    @NotNull
    public Cell createValue(@NotNull String text)
    {
        return createValue(new Paragraph(text));
    }

    @NotNull
    public Cell createValue(@NotNull IBlockElement element)
    {
        Cell cell = createBorderlessBase();
        cell.addStyle(reportResources.valueStyle());
        cell.add(element);
        return cell;
    }

    @NotNull
    private static Cell createBorderlessBase()
    {
        Cell cell = new Cell();
        cell.setBorder(Border.NO_BORDER);
        return cell;
    }
}

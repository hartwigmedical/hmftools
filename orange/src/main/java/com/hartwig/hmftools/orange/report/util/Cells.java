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

public class Cells
{
    private final ReportResources mReportResources;

    public Cells(final ReportResources reportResources)
    {
        mReportResources = reportResources;
    }

    public Cell createHeader(final String text)
    {
        Cell cell = new Cell(1, 1);
        cell.setHeight(23); // Set fixed height to create consistent spacing between table title and header
        cell.setBorder(Border.NO_BORDER);
        cell.setVerticalAlignment(VerticalAlignment.BOTTOM);
        cell.addStyle(mReportResources.tableHeaderStyle());
        cell.add(new Paragraph(text.toUpperCase()));
        return cell;
    }

    public Cell createSpanningEntry(final Table table, final String text)
    {
        Cell cell = new Cell(1, table.getNumberOfColumns());
        cell.add(new Paragraph(text));
        cell.setBorder(Border.NO_BORDER);
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(mReportResources.tableContentStyle());
        return cell;
    }

    public Cell createSpanningWarning(final Table table, final String text)
    {
        Cell cell = new Cell(1, table.getNumberOfColumns());
        cell.add(new Paragraph(text));
        cell.setBorder(Border.NO_BORDER);
        cell.addStyle(mReportResources.qcWarningStyle());
        return cell;
    }

    public Cell createUrl(final String text, final String url)
    {
        Cell cell = createBorderlessBase();
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(mReportResources.urlStyle());
        cell.add(new Paragraph(text));
        cell.setAction(PdfAction.createURI(url));
        return cell;
    }

    public Cell createContent(final String text)
    {
        return createContent(new Paragraph(text));
    }

    public Cell createContent(final IBlockElement element)
    {
        Cell cell = createBorderlessBase();
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(mReportResources.tableContentStyle());
        cell.add(element);
        return cell;
    }

    public Cell createImage(final Image image)
    {
        Cell cell = createBorderlessBase();
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(mReportResources.tableContentStyle());
        cell.add(image);
        return cell;
    }

    public Cell createKey(final String text)
    {
        Cell cell = createBorderlessBase();
        cell.addStyle(mReportResources.keyStyle());
        cell.add(new Paragraph(text));
        return cell;
    }

    public Cell createValue(final String text)
    {
        return createValue(new Paragraph(text));
    }

    public Cell createValue(final IBlockElement element)
    {
        Cell cell = createBorderlessBase();
        cell.addStyle(mReportResources.valueStyle());
        cell.add(element);
        return cell;
    }

    private static Cell createBorderlessBase()
    {
        Cell cell = new Cell();
        cell.setBorder(Border.NO_BORDER);
        return cell;
    }
}

package com.hartwig.hmftools.orange.report.util;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.kernel.pdf.action.PdfAction;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.borders.SolidBorder;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.IBlockElement;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.VerticalAlignment;

import org.jetbrains.annotations.NotNull;

public final class CellUtil {

    private CellUtil() {
    }

    @NotNull
    public static Cell createHeader(@NotNull String text) {
        Cell cell = new Cell(1, 1);
        cell.setHeight(23); // Set fixed height to create consistent spacing between table title and header
        cell.setBorder(Border.NO_BORDER);
        cell.setVerticalAlignment(VerticalAlignment.BOTTOM);
        cell.addStyle(ReportResources.tableHeaderStyle());
        cell.add(new Paragraph(text.toUpperCase()));
        return cell;
    }

    @NotNull
    public static Cell createTransparent(@NotNull String text) {
        Cell cell = createBorderlessBase();
        cell.addStyle(ReportResources.tableContentStyle());
        cell.add(new Paragraph(text));
        return cell;
    }

    @NotNull
    public static Cell createUrl(@NotNull String text, @NotNull String url) {
        Cell cell = createBorderlessBase();
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(ReportResources.urlStyle());
        cell.add(new Paragraph(text));
        cell.setAction(PdfAction.createURI(url));
        return cell;
    }

    @NotNull
    public static Cell createContent(@NotNull String text) {
        return createContent(new Paragraph(text));
    }

    @NotNull
    public static Cell createContent(@NotNull IBlockElement element) {
        Cell cell = createBorderlessBase();
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(ReportResources.tableContentStyle());
        cell.add(element);
        return cell;
    }

    @NotNull
    public static Cell createImage(@NotNull Image image) {
        Cell cell = createBorderlessBase();
        cell.setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25F));
        cell.addStyle(ReportResources.tableContentStyle());
        cell.add(image);
        return cell;
    }

    @NotNull
    public static Cell createKey(@NotNull String text) {
        Cell cell = createBorderlessBase();
        cell.addStyle(ReportResources.keyStyle());
        cell.add(new Paragraph(text));
        return cell;
    }

    @NotNull
    public static Cell createValue(@NotNull String text) {
        Cell cell = createBorderlessBase();
        cell.addStyle(ReportResources.valueStyle());
        cell.add(new Paragraph(text));
        return cell;
    }

    @NotNull
    private static Cell createBorderlessBase() {
        Cell cell = new Cell();
        cell.setBorder(Border.NO_BORDER);
        return cell;
    }
}

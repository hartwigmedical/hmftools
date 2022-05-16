package com.hartwig.hmftools.orange.report.util;

import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Tables {

    private static final float TABLE_BOTTOM_MARGIN = 20;

    private Tables() {
    }

    @NotNull
    public static Table createEmpty(@NotNull String title, float width) {
        Cell headerCell = new Cell().setBorder(Border.NO_BORDER).add(new Paragraph(title).addStyle(ReportResources.tableTitleStyle()));

        Table table = Tables.createContent(width, new float[] { 1 }, new Cell[] { headerCell });
        table.setKeepTogether(true);
        table.setMarginBottom(TABLE_BOTTOM_MARGIN);
        table.addCell(Cells.createContent(new Paragraph(DataUtil.NONE_STRING)));

        return table;
    }

    @NotNull
    public static Table createContent(float width, @NotNull float[] columnPercentageWidths, @NotNull Cell[] headerCells) {
        Table table = new Table(UnitValue.createPercentArray(columnPercentageWidths)).setWidth(width);
        table.setFixedLayout();

        for (Cell headerCell : headerCells) {
            table.addHeaderCell(headerCell);
        }

        return table;
    }

    @NotNull
    public static Table createWrapping(@NotNull Table contentTable) {
        return createWrapping(contentTable, null);
    }

    @NotNull
    public static Table createWrapping(@NotNull Table contentTable, @Nullable String title) {
        contentTable.addFooterCell(new Cell(1, contentTable.getNumberOfColumns()).setBorder(Border.NO_BORDER)
                .setPaddingTop(5)
                .setPaddingBottom(5)
                .add(new Paragraph("The table continues on the next page".toUpperCase()).addStyle(ReportResources.subTextStyle())))
                .setSkipLastFooter(true);

        Table continuedWrapTable = new Table(1).setMinWidth(contentTable.getWidth())
                .addHeaderCell(new Cell().setBorder(Border.NO_BORDER)
                        .add(new Paragraph("Continued from the previous page".toUpperCase()).addStyle(ReportResources.subTextStyle())))
                .setSkipFirstHeader(true)
                .addCell(new Cell().add(contentTable).setPadding(0).setBorder(Border.NO_BORDER));

        Table table = new Table(1).setMinWidth(contentTable.getWidth()).setMarginBottom(TABLE_BOTTOM_MARGIN);
        if (title != null) {
            table.addHeaderCell(new Cell().setBorder(Border.NO_BORDER)
                    .setPadding(0)
                    .add(new Paragraph(title).addStyle(ReportResources.tableTitleStyle())));
        }

        table.addCell(new Cell().add(continuedWrapTable).setPadding(0).setBorder(Border.NO_BORDER));
        return table;
    }
}

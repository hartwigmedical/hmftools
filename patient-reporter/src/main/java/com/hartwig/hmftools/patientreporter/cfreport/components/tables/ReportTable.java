package com.hartwig.hmftools.patientreporter.cfreport.components.tables;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.element.IBlockElement;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import org.jetbrains.annotations.NotNull;

public abstract class ReportTable {

    protected final Table table;

    public ReportTable(float[] columnPercentageWidths, @NotNull String[] headerNames) {

        // Initialize table
        table = new Table(UnitValue.createPercentArray(columnPercentageWidths));
        table.setWidth(ReportResources.CONTENT_WIDTH_WIDE);

        // Add header columns
        for (String name: headerNames) {
            table.addHeaderCell(new TableHeaderCell().add(new Paragraph(name.toUpperCase())));
        }

    }

    public void addCell(String value) {
        addCell(new Paragraph(value));
    }

    public void addCell(IBlockElement element) {
        table.addCell(new TableCell().add(element));
    }

    public Table getTable() {
        return table;
    }

}

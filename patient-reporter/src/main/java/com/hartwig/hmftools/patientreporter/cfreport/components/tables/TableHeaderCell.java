package com.hartwig.hmftools.patientreporter.cfreport.components.tables;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.property.VerticalAlignment;

/**
 * @deprecated
 */
public class TableHeaderCell extends Cell {

    public TableHeaderCell() {
        super();
        setBorder(Border.NO_BORDER);
        setVerticalAlignment(VerticalAlignment.BOTTOM);
        addStyle(ReportResources.tableHeaderStyle());
    }

}

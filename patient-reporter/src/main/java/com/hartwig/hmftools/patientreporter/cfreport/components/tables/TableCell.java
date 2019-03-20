package com.hartwig.hmftools.patientreporter.cfreport.components.tables;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.borders.SolidBorder;
import com.itextpdf.layout.element.Cell;
import org.jetbrains.annotations.NotNull;

/**
 * @deprecated
 */
public class TableCell extends Cell {

    public TableCell() {
        super();
        setBorder(Border.NO_BORDER);
        setBorderBottom(new SolidBorder(ReportResources.PALETTE_MID_GREY, 0.25f));
        addStyle(ReportResources.tableContentStyle());
        setKeepTogether(true);
    }

}

package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;
import org.jetbrains.annotations.NotNull;

public class TumorLocationAndTypeTable {

    public static Table createTumorLocationAndType(@NotNull String primaryTumorLocation, @NotNull String cancerSubType, float width) {

        // Initialize table
        Table table = new Table(UnitValue.createPercentArray(new float[]{1, 1}));
        table.setWidth(width);

        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("PRIMARY TUMOR LOCATION")
                        .addStyle(ReportResources.subTextStyle())));
        table.addCell(TableUtil.getLayoutCell()
                .add(new Paragraph("CANCER SUBTYPE")
                        .addStyle(ReportResources.subTextStyle())));
        table.addCell(TableUtil.getLayoutCell()
                .add(DataLabel.createDataLabel(primaryTumorLocation)));
        table.addCell(TableUtil.getLayoutCell()
                .add(DataLabel.createDataLabel(cancerSubType)));

        return table;

    }

}

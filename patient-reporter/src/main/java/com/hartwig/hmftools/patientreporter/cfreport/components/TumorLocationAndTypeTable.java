package com.hartwig.hmftools.patientreporter.cfreport.components;

import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.UnitValue;

import org.jetbrains.annotations.NotNull;

public final class TumorLocationAndTypeTable {

    private TumorLocationAndTypeTable() {
    }

    @NotNull
    public static Table createBiopsyLocationAndTumorLocation(@NotNull String primaryTumorLocation, @NotNull String biopsyLocation,
            float width) {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 2, 2 }));
        table.setWidth(width);

        table.addCell(TableUtil.createLayoutCell().add(new Paragraph("PRIMARY TUMOR LOCATION").addStyle(ReportResources.subTextStyle())));
        table.addCell(TableUtil.createLayoutCell().add(new Paragraph("BIOPSY LOCATION").addStyle(ReportResources.subTextStyle())));

        table.addCell(TableUtil.createLayoutCell().add(DataLabel.createDataLabel(primaryTumorLocation)));
        table.addCell(TableUtil.createLayoutCell().add(DataLabel.createDataLabel(biopsyLocation)));

        return table;
    }

    @NotNull
    public static Table createTumorType(String primaryTumorType, float width) {
        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 }));
        table.setWidth(width);

        table.addCell(TableUtil.createLayoutCell().add(new Paragraph("PRIMARY TUMOR TYPE").addStyle(ReportResources.subTextStyle())));
        table.addCell(TableUtil.createLayoutCell().add(DataLabel.createDataLabel(primaryTumorType)));

        return table;
    }
}

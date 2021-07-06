package com.hartwig.hmftools.orange.report.components;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.borders.SolidBorder;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Text;
import com.itextpdf.layout.property.BorderRadius;

import org.jetbrains.annotations.NotNull;

public final class DataLabel {

    private DataLabel() {
    }

    @NotNull
    public static Paragraph createDataLabel(@NotNull String text) {
        return new Paragraph().add(new Text(text).addStyle(ReportResources.dataHighlightStyle())
                .setBackgroundColor(ReportResources.PALETTE_BLUE)
                .setFontColor(ReportResources.PALETTE_WHITE)
                .setBorder(new SolidBorder(ReportResources.PALETTE_BLUE, 2))
                .setBorderRadius(new BorderRadius(3f)));
    }
}

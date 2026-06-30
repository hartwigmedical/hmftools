package com.hartwig.hmftools.orange.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import net.sf.dynamicreports.report.builder.style.StyleBuilder;

public final class OrangeStyles
{
    public static final StyleBuilder BORDERED_BOX_STYLE = stl.style()
            .setBorder(stl.border(stl.penThin().setLineColor(OrangeColors.PALETTE_LIGHT_GREY)))
            .setPadding(stl.padding(5));

    public static final StyleBuilder ORANGE_BOX_STYLE = stl.style()
            .setBackgroundColor(OrangeColors.PALETTE_ORANGE)
            .setPadding(stl.padding().setTop(10).setLeft(15).setRight(15));

    private OrangeStyles() {}
}

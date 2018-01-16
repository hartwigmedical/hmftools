package com.hartwig.hmftools.patientreporter.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import java.awt.Color;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.HorizontalTextAlignment;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;

public final class Commons {

    private static final String FONT = "Tinos";
    private static final String MONOSPACE_FONT = "Inconsolata";
    private static final Color BORKIE_COLOR = new Color(221, 235, 247);
    public static final String DATE_TIME_FORMAT = "dd-MMM-yyyy";
    private static final int TABLE_PADDING = 1;
    public static final int SECTION_VERTICAL_GAP = 25;
    public static final int HEADER_TO_DETAIL_VERTICAL_GAP = 8;

    @NotNull
    public static StyleBuilder tableHeaderStyle() {
        return fontStyle().bold()
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setFontSize(10)
                .setBorder(stl.pen1Point())
                .setBackgroundColor(BORKIE_COLOR)
                .setPadding(TABLE_PADDING);
    }

    @NotNull
    public static StyleBuilder dataStyle() {
        return fontStyle().setFontSize(8)
                .setHorizontalTextAlignment(HorizontalTextAlignment.CENTER)
                .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
                .setPadding(TABLE_PADDING);
    }

    @NotNull
    public static StyleBuilder dataTableStyle() {
        return dataStyle().setBorder(stl.pen1Point());
    }

    @NotNull
    static StyleBuilder smallDataTableStyle() {
        return dataStyle().setFontSize(7).setBorder(stl.penThin().setLineColor(Color.black));
    }

    @NotNull
    public static StyleBuilder fontStyle() {
        return stl.style().setFontName(FONT);
    }

    @NotNull
    public static JasperReportBuilder baseTable() {
        return report().setColumnStyle(dataStyle()).setColumnTitleStyle(tableHeaderStyle()).highlightDetailEvenRows();
    }

    @NotNull
    public static JasperReportBuilder monospaceBaseTable() {
        return report().setColumnStyle(dataStyle().setFontName(MONOSPACE_FONT))
                .setColumnTitleStyle(tableHeaderStyle())
                .highlightDetailEvenRows();
    }

    @NotNull
    public static JasperReportBuilder smallTable() {
        return report().setColumnStyle(dataStyle().setFontSize(6)).setColumnTitleStyle(tableHeaderStyle()).highlightDetailEvenRows();
    }

    @NotNull
    public static StyleBuilder sectionHeaderStyle() {
        return fontStyle().bold().setFontSize(12).setHorizontalTextAlignment(HorizontalTextAlignment.CENTER);
    }

    @NotNull
    public static StyleBuilder linkStyle() {
        return dataStyle().setForegroundColor(Color.BLUE);
    }
}

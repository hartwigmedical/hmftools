package com.hartwig.hmftools.orange.report;

import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import net.sf.dynamicreports.report.builder.style.FontBuilder;
import net.sf.dynamicreports.report.builder.style.StyleBuilder;
import net.sf.dynamicreports.report.constant.VerticalTextAlignment;

public final class OrangeFonts
{
    public static final String FONT_FAMILY = "Nimbus Sans";

    private static final int FONT_SIZE_TITLE = 11;
    private static final int FONT_SIZE_CHAPTER_TITLE = 10;
    private static final int FONT_SIZE_TABLE_TITLE = 8;
    private static final int FONT_SIZE_TABLE_HEADER = 7;
    private static final int FONT_SIZE_TABLE_CONTENT = 7;
    private static final int FONT_SIZE_SUB_TEXT = 6;
    private static final int FONT_SIZE_PAGE_NUMBER = 7;
    private static final int FONT_SIZE_SIDE_PANEL_LABEL = 7;
    private static final int FONT_SIZE_SIDE_PANEL_VALUE = 10;

    public static FontBuilder titleFont()
    {
        return stl.font().setFontName(FONT_FAMILY).setBold(true).setFontSize(FONT_SIZE_TITLE);
    }

    public static final StyleBuilder CHAPTER_TITLE_STYLE = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setBold(true).setFontSize(FONT_SIZE_CHAPTER_TITLE))
            .setForegroundColor(OrangeColors.PALETTE_ORANGE);

    public static final StyleBuilder TABLE_TITLE_STYLE = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setBold(true).setFontSize(FONT_SIZE_TABLE_TITLE))
            .setForegroundColor(OrangeColors.PALETTE_ORANGE);

    public static final StyleBuilder TABLE_TITLE_STYLE_WITH_GAP = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setBold(true).setFontSize(FONT_SIZE_TABLE_TITLE))
            .setForegroundColor(OrangeColors.PALETTE_ORANGE)
            .setPadding(stl.padding().setTop(10).setBottom(12));

    public static final StyleBuilder TABLE_HEADER_STYLE_UNPADDED = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setFontSize(FONT_SIZE_TABLE_HEADER))
            .setPadding(stl.padding().setTop(3).setBottom(3).setLeft(0).setRight(0))
            .setForegroundColor(OrangeColors.PALETTE_MID_GREY);

    public static final StyleBuilder TABLE_CONTENT_STYLE_UNPADDED = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setFontSize(FONT_SIZE_TABLE_CONTENT))
            .setPadding(stl.padding().setTop(3).setBottom(3).setLeft(0).setRight(0))
            .setForegroundColor(OrangeColors.PALETTE_DARK_GREY);

    public static final StyleBuilder TABLE_CONTENT_STYLE = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setFontSize(FONT_SIZE_TABLE_CONTENT))
            .setPadding(stl.padding(3))
            .setForegroundColor(OrangeColors.PALETTE_DARK_GREY);

    public static final StyleBuilder SUB_TEXT_STYLE = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setFontSize(FONT_SIZE_SUB_TEXT))
            .setForegroundColor(OrangeColors.PALETTE_BLACK);

    public static final StyleBuilder PAGE_NUMBER_STYLE = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setBold(true).setFontSize(FONT_SIZE_PAGE_NUMBER))
            .setForegroundColor(OrangeColors.PALETTE_ORANGE);

    public static final StyleBuilder SIDE_PANEL_LABEL_STYLE = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setBold(true).setFontSize(FONT_SIZE_SIDE_PANEL_LABEL))
            .setVerticalTextAlignment(VerticalTextAlignment.MIDDLE)
            .setForegroundColor(OrangeColors.PALETTE_WHITE);

    public static final StyleBuilder SIDE_PANEL_VALUE_STYLE = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setBold(true).setFontSize(FONT_SIZE_SIDE_PANEL_VALUE))
            .setForegroundColor(OrangeColors.PALETTE_WHITE);

    public static final StyleBuilder QC_WARNING_STYLE = stl.style()
            .setFont(stl.font().setFontName(FONT_FAMILY).setBold(true).setFontSize(FONT_SIZE_TABLE_CONTENT))
            .setForegroundColor(OrangeColors.PALETTE_DARK_GREY);

    private OrangeFonts()
    {
    }
}

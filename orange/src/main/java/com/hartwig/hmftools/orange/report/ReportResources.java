package com.hartwig.hmftools.orange.report;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import com.hartwig.hmftools.orange.OrangeApplication;
import com.itextpdf.io.font.FontProgram;
import com.itextpdf.io.font.FontProgramFactory;
import com.itextpdf.io.font.PdfEncodings;
import com.itextpdf.kernel.colors.DeviceRgb;
import com.itextpdf.kernel.font.PdfFont;
import com.itextpdf.kernel.font.PdfFontFactory;
import com.itextpdf.layout.Style;

import org.jetbrains.annotations.NotNull;

public final class ReportResources {

    static final String METADATA_TITLE = "HMF ORANGE Report v" + OrangeApplication.VERSION;
    static final String METADATA_AUTHOR = "Hartwig Platinum Pipeline";

    public static final String NOT_AVAILABLE = "NA";

    public static final float PAGE_MARGIN_TOP = 100; // Top margin also excludes the chapter title, which is rendered in the header
    public static final float PAGE_MARGIN_LEFT = 30;
    public static final float PAGE_MARGIN_RIGHT = 30;
    public static final float PAGE_MARGIN_BOTTOM = 40;

    public static final DeviceRgb PALETTE_WHITE = new DeviceRgb(255, 255, 255);
    public static final DeviceRgb PALETTE_BLACK = new DeviceRgb(0, 0, 0);

    public static final DeviceRgb PALETTE_DARK_GREY = new DeviceRgb(39, 47, 50);
    public static final DeviceRgb PALETTE_MID_GREY = new DeviceRgb(101, 106, 108);
    public static final DeviceRgb PALETTE_BLUE = new DeviceRgb(38, 90, 166);

    public static final DeviceRgb PALETTE_ORANGE = new DeviceRgb(242, 139, 31);
    public static final DeviceRgb PALETTE_ORANGE_1 = new DeviceRgb(255, 165, 0);
    public static final DeviceRgb PALETTE_ORANGE_2 = new DeviceRgb(235, 155, 0);
    public static final DeviceRgb PALETTE_ORANGE_3 = new DeviceRgb(215, 145, 0);
    public static final DeviceRgb PALETTE_ORANGE_4 = new DeviceRgb(195, 135, 0);
    public static final DeviceRgb PALETTE_ORANGE_5 = new DeviceRgb(175, 125, 0);
    public static final DeviceRgb PALETTE_ORANGE_6 = new DeviceRgb(155, 115, 0);

    private static final String FONT_REGULAR_PATH = "fonts/nimbus-sans/NimbusSansL-Regular.ttf";
    private static final String FONT_BOLD_PATH = "fonts/nimbus-sans/NimbusSansL-Bold.ttf";

    @NotNull
    public static String formatSingleDigitDecimal(double num) {
        return formatDecimal(num, "0.0");
    }

    @NotNull
    public static String formatTwoDigitDecimal(double num) {
        return formatDecimal(num, "0.00");
    }

    @NotNull
    public static String formatPercentage(double num) {
        return formatPercentage(num, true);
    }

    @NotNull
    public static String formatPercentage(double num, boolean multiplyBy100) {
        return formatDecimal(multiplyBy100 ? num * 100 : num, "0'%'");
    }

    @NotNull
    public static String formatPercentageOneDecimal(double num) {
        return formatDecimal(num * 100, "0.0'%'");
    }

    @NotNull
    private static String formatDecimal(double num, @NotNull String format) {
        // To make sure every decimal format uses a dot as separator rather than a comma.
        return new DecimalFormat(format, DecimalFormatSymbols.getInstance(Locale.ENGLISH)).format(num);
    }

    @NotNull
    public static PdfFont fontRegular() {
        // Cannot be created statically as every PDF needs their own private font objects.
        return createFontFromProgram(loadFontProgram(FONT_REGULAR_PATH));
    }

    @NotNull
    public static PdfFont fontBold() {
        // Cannot be created statically as every PDF needs their own private font objects.
        return createFontFromProgram(loadFontProgram(FONT_BOLD_PATH));
    }

    public static Style chapterTitleStyle() {
        return new Style().setFont(fontBold()).setFontSize(10).setFontColor(ReportResources.PALETTE_ORANGE);
    }

    public static Style tableTitleStyle() {
        return new Style().setFont(fontBold()).setFontSize(8).setFontColor(ReportResources.PALETTE_ORANGE);
    }

    public static Style tableHeaderStyle() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_MID_GREY);
    }

    public static Style tableContentStyle() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_DARK_GREY);
    }

    public static Style keyStyle() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_MID_GREY);
    }

    public static Style valueStyle() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_MID_GREY);
    }

    public static Style subTextStyle() {
        return new Style().setFont(fontRegular()).setFontSize(6).setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style pageNumberStyle() {
        return new Style().setFont(fontBold()).setFontSize(7).setFontColor(ReportResources.PALETTE_ORANGE);
    }

    public static Style sidePanelLabelStyle() {
        return new Style().setFont(fontBold()).setFontSize(7).setFontColor(ReportResources.PALETTE_WHITE);
    }

    public static Style sidePanelValueStyle() {
        return new Style().setFont(fontBold()).setFontSize(10).setFontColor(ReportResources.PALETTE_WHITE);
    }

    public static Style urlStyle() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_BLUE);
    }

    @NotNull
    private static PdfFont createFontFromProgram(@NotNull FontProgram program) {
        return PdfFontFactory.createFont(program, PdfEncodings.IDENTITY_H);
    }

    @NotNull
    private static FontProgram loadFontProgram(@NotNull String resourcePath) {
        try {
            return FontProgramFactory.createFont(resourcePath);
        } catch (IOException exception) {
            // Should never happen, fonts are loaded from code
            throw new IllegalStateException(exception);
        }
    }
}

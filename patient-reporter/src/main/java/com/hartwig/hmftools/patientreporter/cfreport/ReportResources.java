package com.hartwig.hmftools.patientreporter.cfreport;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.time.LocalDate;
import java.util.Locale;

import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.itextpdf.io.font.FontProgram;
import com.itextpdf.io.font.FontProgramFactory;
import com.itextpdf.io.font.PdfEncodings;
import com.itextpdf.kernel.colors.DeviceRgb;
import com.itextpdf.kernel.font.PdfFont;
import com.itextpdf.kernel.font.PdfFontFactory;
import com.itextpdf.layout.Style;

import org.jetbrains.annotations.NotNull;

public final class ReportResources {

    private static final String HARTWIG_NAME = "Hartwig Medical Foundation";
    public static final String HARTWIG_ADDRESS = HARTWIG_NAME + ", Science Park 408, 1098XH Amsterdam";
    public static final String CONTACT_EMAIL_GENERAL = "diagnosticssupport@hartwigmedicalfoundation.nl";
    public static final String CONTACT_EMAIL_QA = "qualitysystem@hartwigmedicalfoundation.nl";
    public static final String SIGNATURE_NAME = "Edwin Cuppen";
    public static final String SIGNATURE_TITLE = "Director " + HARTWIG_NAME;
    public static final String VERSION_REPORT = "version " + PatientReporterApplication.VERSION;
    public static final String MANUAL = "https://www.oncoact.nl/manual";

    public static final double PURITY_CUTOFF = 0.195;

    static final String METADATA_TITLE = "HMF Sequencing Report v" + PatientReporterApplication.VERSION;
    static final String METADATA_AUTHOR = HARTWIG_NAME;
    public static final String REPORT_DATE = DataUtil.formatDate(LocalDate.now());

    static final float PAGE_MARGIN_TOP = 150; // Top margin also excludes the chapter title, which is rendered in the header
    public static final float PAGE_MARGIN_LEFT = 55.5f;
    static final float PAGE_MARGIN_RIGHT = 29;
    static final float PAGE_MARGIN_BOTTOM = 62;

    public static final float CONTENT_WIDTH_NARROW = 330; // Width of the content on a narrow page (page with full side panel)
    public static final float CONTENT_WIDTH_WIDE = 510; // Width of the content on a narrow page (page without full side panel)
    public static final float CONTENT_WIDTH_WIDE_SMALL = 240; // Width of the content on a narrow page (page with full side panel)

    public static final DeviceRgb PALETTE_WHITE = new DeviceRgb(255, 255, 255);
    public static final DeviceRgb PALETTE_BLACK = new DeviceRgb(0, 0, 0);
    public static final DeviceRgb PALETTE_BLUE = new DeviceRgb(38, 90, 166);
    public static final DeviceRgb PALETTE_MID_BLUE = new DeviceRgb(110, 139, 189);
    public static final DeviceRgb PALETTE_DARK_BLUE = new DeviceRgb(93, 85, 164);
    public static final DeviceRgb PALETTE_RED = new DeviceRgb(232, 60, 55);
    public static final DeviceRgb PALETTE_CYAN = new DeviceRgb(0, 179, 233);
    private static final DeviceRgb PALETTE_DARK_GREY = new DeviceRgb(39, 47, 50);
    public static final DeviceRgb PALETTE_MID_GREY = new DeviceRgb(101, 106, 108);
    public static final DeviceRgb PALETTE_LIGHT_GREY = new DeviceRgb(205, 206, 207);
    public static final DeviceRgb PALETTE_PINK = new DeviceRgb(230, 21, 124);
    public static final DeviceRgb PALETTE_VIOLET = new DeviceRgb(156, 97, 168);

    private static final String FONT_REGULAR_PATH = "fonts/nimbus-sans/NimbusSansL-Regular.ttf";
    private static final String FONT_BOLD_PATH = "fonts/nimbus-sans/NimbusSansL-Bold.ttf";
    private static final String ICON_FONT_PATH = "fonts/hmf-icons/hmf-icons.ttf";

    public static final float BODY_TEXT_LEADING = 10F;

    public static float maxPointSizeForWidth(@NotNull PdfFont font, float initialFontSize, float minFontSize, @NotNull String text,
            float maxWidth) {
        float fontIncrement = 0.1f;

        float fontSize = initialFontSize;
        float width = font.getWidth(text, initialFontSize);
        while (width > maxWidth && fontSize > minFontSize) {
            fontSize -= fontIncrement;
            width = font.getWidth(text, fontSize);
        }

        return fontSize;
    }

    @NotNull
    public static DecimalFormat decimalFormat(@NotNull String format) {
        // To make sure every decimal format uses a dot as separator rather than a comma.
        return new DecimalFormat(format, DecimalFormatSymbols.getInstance(Locale.ENGLISH));
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

    @NotNull
    public static PdfFont iconFont() {
        // Cannot be created statically as every PDF needs their own private font objects.
        return createFontFromProgram(loadFontProgram(ICON_FONT_PATH));
    }

    public static Style responseStyle() {
        return new Style().setFont(fontBold()).setFontSize(8).setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style resistentStyle() {
        return new Style().setFont(fontBold()).setFontSize(8).setFontColor(ReportResources.PALETTE_RED);
    }

    public static Style predictedStyle() {
        return new Style().setFont(fontBold()).setFontSize(8).setFontColor(ReportResources.PALETTE_VIOLET);
    }

    public static Style chapterTitleStyle() {
        return new Style().setFont(fontBold()).setFontSize(16).setFontColor(ReportResources.PALETTE_BLUE).setMarginTop(0);
    }

    public static Style sectionTitleStyle() {
        return new Style().setFont(fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style tableHeaderStyle() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_MID_GREY);
    }

    public static Style tableContentStyle() {
        return new Style().setFont(fontRegular()).setFontSize(8).setFontColor(ReportResources.PALETTE_DARK_GREY);
    }

    public static Style bodyTextStyle() {
        return new Style().setFont(fontRegular()).setFontSize(8).setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style smallBodyHeadingStyle() {
        return new Style().setFont(fontBold()).setFontSize(10).setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style smallBodyTextStyle() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style smallBodyTextStyleRed() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_RED);
    }

    public static Style smallBodyBoldTextStyle() {
        return new Style().setFont(fontBold()).setFontSize(7).setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style subTextStyle() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style subTextSmallStyle() {
        return new Style().setFont(fontRegular()).setFontSize(5).setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style subTextBoldStyle() {
        return new Style().setFont(fontBold()).setFontSize(7).setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style dataHighlightStyle() {
        return new Style().setFont(fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style dataHighlightNaStyle() {
        return new Style().setFont(fontBold()).setFontSize(7).setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style pageNumberStyle() {
        return new Style().setFont(fontBold()).setFontSize(8).setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style sidePanelLabelStyle() {
        return new Style().setFont(fontBold()).setFontSize(7).setFontColor(ReportResources.PALETTE_WHITE);
    }

    public static Style sidePanelValueStyle() {
        return new Style().setFont(fontBold()).setFontSize(11).setFontColor(ReportResources.PALETTE_WHITE);
    }

    public static Style dataHighlightLinksStyle() {
        return new Style().setFont(fontRegular()).setFontSize(8).setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style urlStyle() {
        return new Style().setFont(fontRegular()).setFontSize(7).setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style tableTitleStyle() {
        return new Style().setFont(fontBold()).setFontSize(8).setFontColor(ReportResources.PALETTE_BLUE);
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

package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.hartwig.hmftools.patientreporter.cfreport.data.DataUtil;
import com.itextpdf.io.font.FontProgram;
import com.itextpdf.io.font.FontProgramFactory;
import com.itextpdf.io.font.PdfEncodings;
import com.itextpdf.io.image.ImageData;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.colors.DeviceRgb;
import com.itextpdf.kernel.font.PdfFont;
import com.itextpdf.kernel.font.PdfFontFactory;
import com.itextpdf.layout.Style;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.io.IOException;
import java.io.InputStream;
import java.time.LocalDate;

/**
 * Shared info, settings and resources for the report
 */
public final class ReportResources {

    private static final Logger LOGGER = LogManager.getLogger(CFReportWriter.class);

    // PDF Content
    private static final String HARTWIG_NAME = "Hartwig Medical Foundation";
    public static final String HARTWIG_ADDRESS = HARTWIG_NAME + ", Science Park 408, 1098XH Amsterdam";
    public static final String CONTACT_EMAIL_GENERAL = "info@hartwigmedicalfoundation.nl";
    public static final String CONTACT_EMAIL_QA = "qualitysystem@hartwigmedicalfoundation.nl";
    public static final String SIGNATURE_NAME = "Edwin Cuppen";
    public static final String SIGNATURE_TITLE = "Director " + HARTWIG_NAME;

    // PDF Document metadata
    public static final String METADATA_TITLE = "HMF Sequencing Report v" + PatientReporterApplication.VERSION;
    public static final String METADATA_AUTHOR = HARTWIG_NAME;
    public static final String DATE_TIME_FORMAT = "dd-MMM-yyyy";
    public static final String REPORT_DATE = DataUtil.formatDate(LocalDate.now());

    // Page margins for normal content (so excluding header and footer) in pt
    public static final float PAGE_MARGIN_TOP = 150; // Top margin also excludes the chapter title, which is rendered in the header
    public static final float PAGE_MARGIN_LEFT = 55.5f;
    public static final float PAGE_MARGIN_RIGHT = 29;
    public static final float PAGE_MARGIN_BOTTOM = 62;

    public static final float CONTENT_WIDTH_NARROW = 330; // Width of the content on a narrow page (page with full sidepanel)
    public static final float CONTENT_WIDTH_WIDE = 510; // Width of the content on a narrow page (page without full sidepanel)

    // Color palette
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

    // Fonts and text spacing
    private static String FONT_REGULAR_PATH = "fonts/nimbus-sans/NimbusSansL-Regular.ttf";
    private static String FONT_BOLD_PATH = "fonts/nimbus-sans/NimbusSansL-Bold.ttf";
    private static String ICON_FONT_PATH = "fonts/hmf-icons/hmf-icons.ttf";

    private static FontProgram fontProgramRegular = null;
    private static FontProgram fontProgramBold = null;
    private static FontProgram iconFontProgram = null;

    public static final float BODY_TEXT_LEADING = 10f;

    @Nullable
    public static PdfFont getFontRegular() {
        if (fontProgramRegular == null) {
            fontProgramRegular = loadFontProgram(FONT_REGULAR_PATH);
        }
        return createFontFromProgram(fontProgramRegular);
    }

    @Nullable
    public static PdfFont getFontBold() {
        if (fontProgramBold == null) {
            fontProgramBold = loadFontProgram(FONT_BOLD_PATH);
        }
        return createFontFromProgram(fontProgramBold);
    }

    @Nullable
    public static PdfFont getIconFont() {
        if (iconFontProgram == null) {
            iconFontProgram = loadFontProgram(ICON_FONT_PATH);
        }
        return createFontFromProgram(iconFontProgram);
    }

    /**
     * Get the max fitting font size
     * @param font
     * @param initalFontSize
     * @param text
     * @param maxWidth
     * @return smallest font size where text fits in maxWidth, or minFontSize
     */
    public static float getMaxPointSizeForWidth(@NotNull PdfFont font, float initalFontSize, float minFontSize, @NotNull String text, float maxWidth) {
        final float fontIncrement = 0.1f;

        float fontSize = initalFontSize;
        float width = font.getWidth(text, initalFontSize);
        while (width > maxWidth && fontSize > minFontSize) {

            fontSize -= fontIncrement;
            width = font.getWidth(text, fontSize);

        }

        return fontSize;

    }

    public static Style chapterTitleStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(16)
                .setFontColor(ReportResources.PALETTE_BLUE)
                .setMarginTop(0);
    }

    public static Style sectionTitleStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(11)
                .setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style tableHeaderStyle() {
        return new Style()
                .setFont(getFontRegular())
                .setFontSize(7)
                .setFontColor(ReportResources.PALETTE_MID_GREY);

    }

    public static Style tableContentStyle() {
        return new Style()
                .setFont(getFontRegular())
                .setFontSize(8)
                .setFontColor(ReportResources.PALETTE_DARK_GREY);
    }

    public static Style bodyTextStyle() {
        return new Style()
                .setFont(getFontRegular())
                .setFontSize(8)
                .setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style smallBodyHeadingStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(10)
                .setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style smallBodyTextStyle() {
        return new Style()
                .setFont(getFontRegular())
                .setFontSize(7)
                .setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style smallBodyBoldTextStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(7)
                .setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style subTextStyle() {
        return new Style()
                .setFont(getFontRegular())
                .setFontSize(7)
                .setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style subTextBoldStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(7)
                .setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static Style dataHighlightStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(11)
                .setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style dataHighlightNaStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(7)
                .setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style pageNumberStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(8)
                .setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static Style sidepanelLabelStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(7)
                .setFontColor(ReportResources.PALETTE_WHITE);
    }

    public static Style sidepanelValueStyle() {
        return new Style()
                .setFont(getFontBold())
                .setFontSize(11)
                .setFontColor(ReportResources.PALETTE_WHITE);
    }

    /**
     * Load image data from resource stream. Returns null when loading fails
     *
     * @param resourcePath
     */
    public static ImageData loadImageData(@NotNull String resourcePath) {

        try {

            byte[] data = loadResourceData(resourcePath);
            if (data == null) {
                throw new Exception("Failed to load image data from " + resourcePath);
            }

            return ImageDataFactory.create(data, true);

        } catch (Exception e) {

            LOGGER.warn(e.getMessage());
            return null;

        }

    }

    /**
     * Load byte array from resource
     *
     * @param resourcePath
     * @throws IOException
     */
    private static byte[] loadResourceData(String resourcePath) throws IOException {

        byte[] data = null;
        InputStream is = ClassLoader.getSystemResourceAsStream(resourcePath);
        if (is != null) {
            data = new byte[is.available()];
            is.read(data);
        }

        return data;

    }

    @Nullable
    private static FontProgram loadFontProgram(@NotNull String resourcePath) {
        try {
            return FontProgramFactory.createFont(resourcePath);
        } catch (Exception e) {
            LOGGER.warn(e.getMessage());
            return null;
        }
    }

    private static PdfFont createFontFromProgram(@NotNull FontProgram program) {
        return PdfFontFactory.createFont(program, PdfEncodings.IDENTITY_H);
    }

}

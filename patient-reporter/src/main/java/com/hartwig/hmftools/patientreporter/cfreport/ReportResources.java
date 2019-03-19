package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.itextpdf.io.image.ImageData;
import com.itextpdf.io.image.ImageDataFactory;
import com.itextpdf.kernel.colors.DeviceRgb;
import com.itextpdf.kernel.font.PdfFont;
import com.itextpdf.layout.Style;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;
import java.io.InputStream;

/**
 * Shared info, settings and resources for the report
 */
public final class ReportResources {

    private static final Logger LOGGER = LogManager.getLogger(CFReportWriter.class);

    public static final String HARTWIG_NAME = "Hartwig Medical Foundation";
    public static final String HARTWIG_ADDRESS = HARTWIG_NAME + ", Science Park 408, 1098XH Amsterdam";

    // PDF Document metadata
    public static final String METADATA_TITLE = "HMF Sequencing Report v" + PatientReporterApplication.VERSION;
    public static final String METADATA_AUTHOR = HARTWIG_NAME;

    // Page margins for normal content (so excluding header and footer) in pt
    public static final float PAGE_MARGIN_TOP = 115;
    public static final float PAGE_MARGIN_LEFT = 75;
    public static final float PAGE_MARGIN_RIGHT = 45;
    public static final float PAGE_MARGIN_BOTTOM = 100;

    // Color palette
    public static final DeviceRgb PALETTE_BLUE = new DeviceRgb(38, 90, 166);
    public static final DeviceRgb PALETTE_MID_BLUE = new DeviceRgb(110, 139, 189);
    public static final DeviceRgb PALETTE_RED = new DeviceRgb(232, 60, 55);
    public static final DeviceRgb PALETTE_CYAN = new DeviceRgb(0, 179, 233);
    public static final DeviceRgb PALETTE_DARK_GREY = new DeviceRgb(39, 47, 50);
    public static final DeviceRgb PALETTE_MID_GREY = new DeviceRgb(101, 106, 108);
    public static final DeviceRgb PALETTE_LIGHT_GREY = new DeviceRgb(205, 206, 207);
    public static final DeviceRgb PALETTE_PINK = new DeviceRgb(230, 21, 124);
    public static final DeviceRgb PALETTE_BLACK = new DeviceRgb(0, 0, 0);

    // Fonts
    private static PdfFont mainFontBook = null;
    private static PdfFont mainFontBold = null;

    private ReportResources() {}

    public static final PdfFont getMainFontBook() {
        if (mainFontBook == null) {
            // Load font
        }
        return mainFontBook;
    }

    public static final PdfFont getMainFontBold() {
        if (mainFontBold == null) {
            // Load font
        }
        return mainFontBold;
    }

    public static final Style chapterTitleStyle() {
        return new Style()
                .setFont(getMainFontBold())
                .setFontSize(16)
                .setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static final Style sectionTitleStyle() {
        return new Style()
                .setFont(getMainFontBold())
                .setFontSize(11)
                .setFontColor(ReportResources.PALETTE_BLUE);
    }

    public static final Style tableheaderStyle() {
        return new Style()
                .setFont(getMainFontBook())
                .setFontSize(7)
                .setFontColor(ReportResources.PALETTE_DARK_GREY);

    }

    public static final Style tableContentStyle() {
        return new Style()
                .setFont(getMainFontBook())
                .setFontSize(8)
                .setFontColor(ReportResources.PALETTE_BLACK);
    }

    public static final Style subTextStyle() {
        return new Style()
                .setFont(getMainFontBook())
                .setFontSize(7)
                .setFontColor(ReportResources.PALETTE_BLACK);
    }

    /**
     * Load image data from resource stream. Returns null when loading fails
     *
     * @param resourcePath
     * @return
     */
    public static final ImageData loadImageData(@NotNull String resourcePath) {

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
     * @return
     * @throws IOException
     */
    public static final byte[] loadResourceData(String resourcePath) throws IOException {

        byte[] data = null;
        InputStream is = new ReportResources().getClass().getClassLoader().getResourceAsStream(resourcePath);
        data = new byte[is.available()];
        is.read(data);

        return data;

    }

}

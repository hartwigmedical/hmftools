package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.PatientReporterApplication;
import com.itextpdf.kernel.colors.DeviceRgb;

/**
 * Static configuration variables for the report
 */
public class ReportConfiguration {

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

}

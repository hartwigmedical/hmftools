package com.hartwig.hmftools.patientreporter.cfreport;

import com.hartwig.hmftools.patientreporter.PatientReporterApplication;

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

}

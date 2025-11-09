package com.hartwig.hmftools.qsee.common;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

public class QseeFileCommon
{
    public static final String QSEE_FILE_ID = "qsee";

    public static final String COL_SAMPLE_ID = "SampleId";
    public static final String COL_SAMPLE_TYPE = "SampleType";
    public static final String COL_FEATURE_TYPE = "FeatureType";
    public static final String COL_FEATURE_NAME = "FeatureName";
    public static final String COL_SOURCE_TOOL = "SourceTool";

    public static final DecimalFormat DECIMAL_FORMAT = QseeFileCommon.createDecimalFormat();

    private static DecimalFormat createDecimalFormat()
    {
        DecimalFormatSymbols symbols = new DecimalFormatSymbols();
        symbols.setInfinity("Inf");
        return new DecimalFormat("0.########", symbols);
    }
}

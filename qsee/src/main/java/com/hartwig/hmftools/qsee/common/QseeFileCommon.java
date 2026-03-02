package com.hartwig.hmftools.qsee.common;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;
import java.util.StringJoiner;

import org.jetbrains.annotations.Nullable;

public class QseeFileCommon
{
    public static final String QSEE_FILE_ID = "qsee";

    public static final String COL_SAMPLE_ID = "SampleId";
    public static final String COL_SAMPLE_TYPE = "SampleType";
    public static final String COL_FEATURE_TYPE = "FeatureType";
    public static final String COL_FEATURE_NAME = "FeatureName";
    public static final String COL_SOURCE_TOOL = "SourceTool";
    public static final String COL_FEATURE_VALUE = "FeatureValue";

    public static final DecimalFormat DECIMAL_FORMAT = QseeFileCommon.createDecimalFormat();

    private static DecimalFormat createDecimalFormat()
    {
        DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.ENGLISH);
        symbols.setInfinity("Inf");
        return new DecimalFormat("0.########", symbols);
    }

    public static String generateFilename(String basePath, String sampleId, String fileId, @Nullable String outputId, String fileExtension)
    {
        StringJoiner filename = new StringJoiner(".");
        filename.add(sampleId);
        filename.add(QSEE_FILE_ID);
        filename.add(fileId);

        if(outputId != null)
            filename.add(outputId);

        filename.add(fileExtension);

        return checkAddDirSeparator(basePath) + filename;
    }

    public static String generateCohortFilename(String basePath, String fileId, @Nullable String outputId, String fileExtension)
    {
        StringJoiner filename = new StringJoiner(".");
        filename.add(QSEE_FILE_ID);
        filename.add("cohort");
        filename.add(fileId);

        if(outputId != null)
            filename.add(outputId);

        filename.add(fileExtension);

        return checkAddDirSeparator(basePath) + filename;
    }
}

package com.hartwig.hmftools.common.basequal.jitter;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.jetbrains.annotations.NotNull;

public class JitterModelParamsFile
{
    enum Column
    {
        unit,
        optimalScaleRepeat4,
        optimalScaleRepeat5,
        optimalScaleRepeat6,
        scaleFitGradient,
        scaleFitIntercept,
        microsatelliteSkew,
    }

    private static String FILE_EXTENSION = ".errorprofile.jitter_params.tsv";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(final String filename, @NotNull final List<JitterModelParams> jitterModelParamsList)
    {
        List<String> columns = Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList());

        DelimFileWriter.write(filename, columns, jitterModelParamsList, (jitterParams, row) ->
        {
            row.set(Column.unit, jitterParams.RepeatUnit);
            row.set(Column.optimalScaleRepeat4, jitterParams.OptimalScaleRepeat4);
            row.set(Column.optimalScaleRepeat5, jitterParams.OptimalScaleRepeat5);
            row.set(Column.optimalScaleRepeat6, jitterParams.OptimalScaleRepeat6);
            row.set(Column.scaleFitGradient, jitterParams.ScaleFitGradient);
            row.set(Column.scaleFitIntercept, jitterParams.ScaleFitIntercept);
            row.set(Column.microsatelliteSkew, jitterParams.MicrosatelliteSkew);
        });
    }

}

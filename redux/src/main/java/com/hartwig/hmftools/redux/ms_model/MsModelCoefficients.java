package com.hartwig.hmftools.redux.ms_model;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.redux.ReduxCommon.COL_NUM_UNITS;
import static com.hartwig.hmftools.common.redux.ReduxCommon.COL_UNIT;
import static com.hartwig.hmftools.common.redux.ReduxCommon.REDUX_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public class MsModelCoefficients
{
    public final String RepeatUnit;
    public final int RepeatCount;
    public final double Coefficient1;
    public final double Coefficient2;

    public MsModelCoefficients(final String repeatUnit, final int repeatCount, final double coefficient1, final double coefficient2)
    {
        RepeatUnit = repeatUnit;
        RepeatCount = repeatCount;
        Coefficient1 = coefficient1;
        Coefficient2 = coefficient2;
    }

    public static String generateFilename(final String basePath, @Nullable final String fileId)
    {
        String filename = basePath + File.separator + "ms_model_coefficients";

        if(fileId != null)
            filename += "." + fileId;

        return filename + TSV_EXTENSION;
    }

    private static final String COL_COEF_1 = "coefficient1";
    private static final String COL_COEF_2 = "coefficient2";

    public static void write(final String filename, final List<MsModelCoefficients> modelCoefficients)
    {
        try
        {
            List<String> lines = Lists.newArrayList();

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(COL_UNIT);
            sj.add(COL_NUM_UNITS);
            sj.add(COL_COEF_1);
            sj.add(COL_COEF_2);

            lines.add(sj.toString());

            for(MsModelCoefficients modelCoeffs : modelCoefficients)
            {
                sj = new StringJoiner(TSV_DELIM);
                sj.add(modelCoeffs.RepeatUnit);
                sj.add(String.valueOf(modelCoeffs.RepeatCount));
                sj.add(String.valueOf(modelCoeffs.Coefficient1));
                sj.add(String.valueOf(modelCoeffs.Coefficient2));
                lines.add(sj.toString());
            }

            Files.write(new File(filename).toPath(), lines);
        }
        catch(IOException e)
        {
            RD_LOGGER.error("failed to write MS model coefficients file({}): {}", filename, e.toString());
        }
    }

    public static List<MsModelCoefficients> read(final String filename)
    {
        List<MsModelCoefficients> modelCoefficients = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            String header = lines.get(0);
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                modelCoefficients.add(new MsModelCoefficients(
                        values[fieldsIndexMap.get(COL_UNIT)],
                        Integer.parseInt(values[fieldsIndexMap.get(COL_NUM_UNITS)]),
                        Double.parseDouble(values[fieldsIndexMap.get(COL_COEF_1)]),
                        Double.parseDouble(values[fieldsIndexMap.get(COL_COEF_2)])));
            }

            return modelCoefficients;
        }
        catch(IOException e)
        {
            RD_LOGGER.error("failed to load MS model coefficients file({}): {}", filename, e.toString());
            return null;
        }
    }
}

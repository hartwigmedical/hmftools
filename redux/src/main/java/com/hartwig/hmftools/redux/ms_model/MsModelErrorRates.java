package com.hartwig.hmftools.redux.ms_model;

import static com.hartwig.hmftools.common.redux.ReduxCommon.COL_NUM_UNITS;
import static com.hartwig.hmftools.common.redux.ReduxCommon.COL_UNIT;
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

public class MsModelErrorRates
{
    public final String RepeatUnit;
    public final int RepeatCount;
    public final double ErrorRate;

    public MsModelErrorRates(final String repeatUnit, final int repeatCount, final double errorRate)
    {
        RepeatUnit = repeatUnit;
        RepeatCount = repeatCount;
        ErrorRate = errorRate;
    }

    public static String generateFilename(final String basePath, @Nullable final String fileId)
    {
        String filename = basePath + File.separator + "ms_model_error_rates";

        if(fileId != null)
            filename += "." + fileId;

        return filename + TSV_EXTENSION;
    }

    private static final String COL_ERROR_RATE = "errorRate";

    public static void write(final String filename, final List<MsModelErrorRates> modelCoefficients)
    {
        try
        {
            List<String> lines = Lists.newArrayList();

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(COL_UNIT);
            sj.add(COL_NUM_UNITS);
            sj.add(COL_ERROR_RATE);

            lines.add(sj.toString());

            for(MsModelErrorRates errorRates : modelCoefficients)
            {
                sj = new StringJoiner(TSV_DELIM);
                sj.add(errorRates.RepeatUnit);
                sj.add(String.valueOf(errorRates.RepeatCount));
                sj.add(String.valueOf(errorRates.ErrorRate));
                lines.add(sj.toString());
            }

            Files.write(new File(filename).toPath(), lines);
        }
        catch(IOException e)
        {
            RD_LOGGER.error("failed to write MS model error rates file({}): {}", filename, e.toString());
        }
    }

    public static List<MsModelErrorRates> read(final String filename)
    {
        List<MsModelErrorRates> modelErrorRates = Lists.newArrayList();

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            String header = lines.get(0);
            Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
            lines.remove(0);

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM, -1);

                modelErrorRates.add(new MsModelErrorRates(
                        values[fieldsIndexMap.get(COL_UNIT)],
                        Integer.parseInt(values[fieldsIndexMap.get(COL_NUM_UNITS)]),
                        Double.parseDouble(values[fieldsIndexMap.get(COL_ERROR_RATE)])));
            }

            return modelErrorRates;
        }
        catch(IOException e)
        {
            RD_LOGGER.error("failed to load MS model error rates file({}): {}", filename, e.toString());
            return null;
        }
    }
}

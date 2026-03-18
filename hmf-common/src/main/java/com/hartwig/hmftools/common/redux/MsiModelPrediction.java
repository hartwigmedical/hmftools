package com.hartwig.hmftools.common.redux;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.redux.ReduxCommon.LOGGER;
import static com.hartwig.hmftools.common.redux.ReduxCommon.REDUX_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

public class MsiModelPrediction
{
    public final double PredictedMsiIndelsPerMb;

    public MsiModelPrediction(double predictedMsiIndelsPerMb)
    {
        PredictedMsiIndelsPerMb = predictedMsiIndelsPerMb;
    }

    public static double INVALID_VALUE = -1;

    public static final MsiModelPrediction INVALID_PREDICTION = new MsiModelPrediction(INVALID_VALUE);

    private static String FILE_EXTENSION = REDUX_FILE_ID + ".msi_prediction.tsv";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(final String filename, final double predictedMsiIndelsPerMb) throws IOException
    {
        List<String> lines = Lists.newArrayListWithCapacity(2);
        lines.add(FLD_PREDICTED_VALUE);
        lines.add(format("%.4f", predictedMsiIndelsPerMb));
        Files.write(new File(filename).toPath(), lines);
    }

    private static String FLD_PREDICTED_VALUE = "PredictedValue";

    public static MsiModelPrediction read(final String filename) throws IOException
    {
        try
        {
            List<String> lines = Files.readAllLines(new File(filename).toPath());

            String header = lines.get(0);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            String[] values = lines.get(0).split(TSV_DELIM, -1);

            double predictedValue = Double.parseDouble(values[fieldsIndexMap.get(FLD_PREDICTED_VALUE)]);

            return new MsiModelPrediction(predictedValue);
        }
        catch(Exception e)
        {
            LOGGER.error("failed to read MSI prediction file({}): {}", filename, e.toString());
            return INVALID_PREDICTION;
        }
    }
}

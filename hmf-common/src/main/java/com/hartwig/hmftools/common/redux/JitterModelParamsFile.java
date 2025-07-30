package com.hartwig.hmftools.common.redux;

import static com.hartwig.hmftools.common.redux.ReduxCommon.REDUX_FILE_ID;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sequencing.ConsensusType;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public class JitterModelParamsFile
{
    enum Column
    {
        unit,
        consensusType,
        optimalScaleRepeat4,
        optimalScaleRepeat5,
        optimalScaleRepeat6,
        scaleFitGradient,
        scaleFitIntercept,
        microsatelliteSkew,
    }

    private static String FILE_EXTENSION = REDUX_FILE_ID + ".jitter_params.tsv";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(final String filename, final List<JitterModelParams> jitterModelParamsList)
    {
        List<String> columns = Arrays.stream(Column.values()).map(Enum::name).collect(Collectors.toList());

        DelimFileWriter.write(filename, columns, jitterModelParamsList, (jitterParams, row) ->
        {
            row.set(Column.unit, jitterParams.RepeatUnit);
            row.set(Column.consensusType, jitterParams.ConsensusType.name());
            row.set(Column.optimalScaleRepeat4, jitterParams.OptimalScaleRepeat4);
            row.set(Column.optimalScaleRepeat5, jitterParams.OptimalScaleRepeat5);
            row.set(Column.optimalScaleRepeat6, jitterParams.OptimalScaleRepeat6);
            row.set(Column.scaleFitGradient, jitterParams.ScaleFitGradient);
            row.set(Column.scaleFitIntercept, jitterParams.ScaleFitIntercept);
            row.set(Column.microsatelliteSkew, jitterParams.MicrosatelliteSkew);
        });
    }

    public static List<JitterModelParams> read(final String filename) throws IOException
    {
        List<JitterModelParams> jitterModelParams = Lists.newArrayList();

        DelimFileReader reader = new DelimFileReader(filename);

        boolean hasConsensusType = reader.hasColumn(Column.consensusType);

        for(DelimFileReader.Row row : reader)
        {
            JitterModelParams params = new JitterModelParams(
                    row.get(Column.unit),
                    hasConsensusType ? ConsensusType.valueOf(row.get(Column.consensusType)) : ConsensusType.IGNORE,
                    row.getDouble(Column.optimalScaleRepeat4),
                    row.getDouble(Column.optimalScaleRepeat5),
                    row.getDouble(Column.optimalScaleRepeat6),
                    row.getDouble(Column.scaleFitGradient),
                    row.getDouble(Column.scaleFitIntercept),
                    row.getDouble(Column.microsatelliteSkew));

            jitterModelParams.add(params);
        }

        return jitterModelParams;
    }
}

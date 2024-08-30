package com.hartwig.hmftools.common.vchord;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.commons.lang3.Validate;

public class VChordPredictionFile
{
    enum Column
    {
        sampleId,
        breastCancerHrdScore,
        ovarianCancerHrdScore,
        pancreaticCancerScore,
        prostateCancerScore,
        otherCancerScore
    }


    private static final String FILE_EXTENSION = ".vchord.prediction.tsv";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static VChordPrediction read(final String filename)
    {
        try (DelimFileReader reader = new DelimFileReader(filename))
        {
            // there should only be one row
            List<VChordPrediction> vChordPredictions = reader.stream().map(row -> ImmutableVChordPrediction.builder()
                    .breastCancerHrdScore(row.getDouble(Column.breastCancerHrdScore))
                    .ovarianCancerHrdScore(row.getDouble(Column.ovarianCancerHrdScore))
                    .pancreaticCancerScore(row.getDouble(Column.pancreaticCancerScore))
                    .prostateCancerScore(row.getDouble(Column.prostateCancerScore))
                    .otherCancerScore(row.getDouble(Column.otherCancerScore))
                    .build()).collect(Collectors.toList());

            Validate.isTrue(vChordPredictions.size() == 1, "Number of vchord prediction record must be 1");

            return vChordPredictions.get(0);
        }
    }

    public static void write(final String filename, final String sampleId, VChordPrediction vChordPrediction)
    {
        try (DelimFileWriter<VChordPrediction> writer = new DelimFileWriter<>(filename, Column.values(),
                (tl, row) ->
                {
                    row.set(Column.sampleId, sampleId);
                    row.set(Column.breastCancerHrdScore, vChordPrediction.breastCancerHrdScore());
                    row.set(Column.ovarianCancerHrdScore, vChordPrediction.ovarianCancerHrdScore());
                    row.set(Column.pancreaticCancerScore, vChordPrediction.pancreaticCancerScore());
                    row.set(Column.prostateCancerScore, vChordPrediction.prostateCancerScore());
                    row.set(Column.otherCancerScore, vChordPrediction.otherCancerScore());
                }))
        {
            writer.writeRow(vChordPrediction);
        }
    }
}

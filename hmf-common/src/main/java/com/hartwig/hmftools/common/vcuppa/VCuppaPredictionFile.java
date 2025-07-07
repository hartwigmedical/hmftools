package com.hartwig.hmftools.common.vcuppa;

import java.io.File;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public class VCuppaPredictionFile
{
    enum Column
    {
        cancerType,
        probability
    }


    private static final String FILE_EXTENSION = ".vcuppa.prediction.tsv";

    public static String generateFilename(String basePath, String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    /*
    public static VChordPrediction read(final String filename)
    {
        try(DelimFileReader reader = new DelimFileReader(filename))
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
    }*/

    public static void write(final String filename, final String sampleId, VCuppaPrediction vCuppaPrediction)
    {
        try(DelimFileWriter<VCuppaPrediction.CancerTypePrediction> writer = new DelimFileWriter<>(filename, Column.values(),
                (cancerTypePrediction, row) ->
                {
                    row.set(Column.cancerType, cancerTypePrediction.cancerType());
                    row.set(Column.probability, cancerTypePrediction.probability());
                }))
        {
            vCuppaPrediction.cancerTypePredictions().forEach(writer::writeRow);
        }
    }
}

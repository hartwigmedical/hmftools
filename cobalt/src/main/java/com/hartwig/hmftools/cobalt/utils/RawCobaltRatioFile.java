package com.hartwig.hmftools.cobalt.utils;

import static com.hartwig.hmftools.common.cobalt.CobaltRatioFile.Column.*;
import static com.hartwig.hmftools.common.cobalt.CobaltRatioFile.FORMAT;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public class RawCobaltRatioFile
{
    private final String fileName;

    public static void write(final String fileName, Collection<RawCobaltRatio> ratios) throws IOException
    {
        CobaltRatioFile.Column[] columns = new CobaltRatioFile.Column[7];
        columns[0] = chromosome;
        columns[1] = position;
        columns[2] = referenceReadDepth;
        columns[3] = tumorReadDepth;
        columns[4] = referenceGCRatio;
        columns[5] = tumorGCRatio;
        columns[6] = referenceGCDiploidRatio;
        try(BufferedWriter writer = createGzipBufferedWriter(fileName))
        {
            DelimFileWriter.write(writer, columns, ratios,
                    (ratio, row) -> {
                        row.set(chromosome, ratio.chromosome());
                        row.set(position, ratio.position());
                        row.set(referenceReadDepth, ratio.referenceReadCount(), FORMAT);
                        row.set(tumorReadDepth, ratio.tumorReadCount(), FORMAT);
                        row.set(referenceGCRatio, ratio.referenceGcRatio(), FORMAT);
                        row.set(tumorGCRatio, ratio.tumorGcRatio(), FORMAT);
                        row.set(referenceGCDiploidRatio, ratio.referenceGcDiploidRatio(), FORMAT);
                    });
        }
    }

    public RawCobaltRatioFile(final String fileName)
    {
        this.fileName = fileName;
    }

    List<RawCobaltRatio> read()
    {
        List<RawCobaltRatio> ratios = new ArrayList<>();
        try(DelimFileReader reader = new DelimFileReader(fileName))
        {
            Integer chrIndex = reader.getColumnIndex(chromosome);
            Integer posIndex = reader.getColumnIndex(position);
            Integer refReadCountIndex = reader.getColumnIndex(referenceReadDepth, "referenceReadCount");
            Integer tumorReadCountIndex = reader.getColumnIndex(tumorReadDepth, "tumorReadCount");
            Integer refGcRatioIndex = reader.getColumnIndex(referenceGCRatio);
            Integer tumorGcRatioIndex = reader.getColumnIndex(tumorGCRatio);
            Integer refGcDiploidRatioIndex = reader.getColumnIndex(referenceGCDiploidRatio);

            for(DelimFileReader.Row row : reader)
            {
                RawCobaltRatio ratio = new RawCobaltRatio(
                        row.get(chrIndex),
                        row.getInt(posIndex),
                        row.getDouble(refReadCountIndex),
                        row.getDouble(tumorReadCountIndex),
                        row.getDouble(refGcRatioIndex),
                        row.getDouble(tumorGcRatioIndex),
                        row.getDouble(refGcDiploidRatioIndex));

                ratios.add(ratio);
            }
        }
        return ratios;
    }
}

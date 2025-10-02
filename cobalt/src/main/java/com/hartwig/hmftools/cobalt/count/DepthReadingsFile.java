package com.hartwig.hmftools.cobalt.count;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Locale;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public class DepthReadingsFile
{
    public enum Column
    {
        chromosome,
        position,
        readDepth,
        gcRatio,
    }

    public static final DecimalFormat FORMAT = new DecimalFormat("#.####", new DecimalFormatSymbols(Locale.ENGLISH));

    public static String generateFilename(final String basePath, final String sample, boolean tumor)
    {
        String infix = tumor ? ".tumor" : ".reference";
        return checkAddDirSeparator(basePath) + sample + infix + ".depths.tsv.gz";
    }

    private final String fileName;

    public static void write(final String fileName, Collection<DepthReading> ratios) throws IOException
    {
        DepthReadingsFile.Column[] columns = new DepthReadingsFile.Column[4];
        columns[0] = Column.chromosome;
        columns[1] = Column.position;
        columns[2] = Column.readDepth;
        columns[3] = Column.gcRatio;
        try(BufferedWriter writer = createGzipBufferedWriter(fileName))
        {
            DelimFileWriter.write(writer, columns, ratios,
                    (ratio, row) ->
                    {
                        row.set(Column.chromosome, ratio.Chromosome);
                        row.set(Column.position, ratio.StartPosition);
                        row.set(Column.readDepth, ratio.ReadDepth, FORMAT);
                        row.set(Column.gcRatio, ratio.ReadGcContent, FORMAT);
                    });
        }
    }

    public DepthReadingsFile(final String fileName)
    {
        this.fileName = fileName;
    }

    List<DepthReading> read()
    {
        List<DepthReading> ratios = new ArrayList<>();
        try(DelimFileReader reader = new DelimFileReader(fileName))
        {
            Integer chrIndex = reader.getColumnIndex(Column.chromosome);
            Integer posIndex = reader.getColumnIndex(Column.position);
            Integer refReadCountIndex = reader.getColumnIndex(Column.readDepth);
            Integer tumorReadCountIndex = reader.getColumnIndex(Column.gcRatio);

            for(DelimFileReader.Row row : reader)
            {
                DepthReading ratio = new DepthReading(
                        row.get(chrIndex),
                        row.getInt(posIndex),
                        row.getDouble(refReadCountIndex),
                        row.getDouble(tumorReadCountIndex));
                ratios.add(ratio);
            }
        }
        return ratios;
    }
}

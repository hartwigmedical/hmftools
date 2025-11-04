package com.hartwig.hmftools.cobalt.metrics;

import static com.hartwig.hmftools.common.cobalt.CobaltRatioFile.FORMAT;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

@SuppressWarnings("DataFlowIssue")
public class WindowStatisticsFile
{
    public enum Column
    {
        chromosome,
        position,
        count,
        mean,
        median,
        min,
        max,
        sd
    }

    public static String fileName(String outputDir, String sample)
    {
        return checkAddDirSeparator(outputDir) + sample + ".fragment_length.statistics.tsv.gz";
    }

    public static void write(final String fileName, List<WindowStatistics> statistics) throws IOException
    {
        WindowStatisticsFile.Column[] columns = WindowStatisticsFile.Column.values();
        try(BufferedWriter writer = createGzipBufferedWriter(fileName))
        {
            DelimFileWriter.write(writer, columns, statistics,
                    (windowStatistics, row) ->
                    {
                        row.set(Column.chromosome, windowStatistics.chromosome());
                        row.set(Column.position, windowStatistics.position());
                        row.set(Column.count, windowStatistics.count());
                        row.set(Column.mean, windowStatistics.mean(), FORMAT);
                        row.set(Column.median, windowStatistics.median(), FORMAT);
                        row.set(Column.min, windowStatistics.min(), FORMAT);
                        row.set(Column.max, windowStatistics.max(), FORMAT);
                        row.set(Column.sd, windowStatistics.sd(), FORMAT);
                    });
        }
    }

    public static List<WindowStatistics> read(String fileName)
    {
        List<WindowStatistics> ratios = new ArrayList<>();
        try(DelimFileReader reader = new DelimFileReader(fileName))
        {
            Integer chrIndex = reader.getColumnIndex(Column.chromosome);
            Integer posIndex = reader.getColumnIndex(Column.position);
            Integer meanIndex = reader.getColumnIndex(Column.mean);
            Integer medianIndex = reader.getColumnIndex(Column.median);
            Integer countIndex = reader.getColumnIndex(Column.count);
            Integer minIndex = reader.getColumnIndex(Column.min);
            Integer maxIndex = reader.getColumnIndex(Column.max);
            Integer sdIndex = reader.getColumnIndex(Column.sd);

            for(DelimFileReader.Row row : reader)
            {
                WindowStatistics ratio = new WindowStatistics(
                        row.get(chrIndex),
                        row.getInt(posIndex),
                        row.getInt(countIndex),
                        row.getDouble(meanIndex),
                        row.getDouble(medianIndex),
                        row.getDouble(minIndex),
                        row.getDouble(maxIndex),
                        row.getDouble(sdIndex));

                ratios.add(ratio);
            }
        }
        return ratios;
    }
}

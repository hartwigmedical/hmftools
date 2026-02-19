package com.hartwig.hmftools.amber.blacklist;

import java.io.File;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public class AmberBlacklistFile
{
    private static final DecimalFormat FORMAT = new DecimalFormat("0.00", new DecimalFormatSymbols(Locale.ENGLISH));
    private static final List<String> COLUMNS = List.of("Chromosome", "Position", "Count", "MeanVaf");

    public static void writeToFile(File destination, List<AmberBlacklistPoint> statistics)
    {
        DelimFileWriter.write(destination.getAbsolutePath(), COLUMNS, statistics, (statistic, row) ->
        {
            row.set(COLUMNS.get(0), statistic.chromosome());
            row.set(COLUMNS.get(1), statistic.position());
            row.set(COLUMNS.get(2), statistic.count());
            row.set(COLUMNS.get(3), FORMAT.format(statistic.meanVaf()));
        });
    }

    public static List<AmberBlacklistPoint> readFromFile(File file)
    {
        try(DelimFileReader reader = new DelimFileReader(file.getAbsolutePath()))
        {
            return reader.stream()
                    .map(row -> new AmberBlacklistPoint(row.get(0), row.getInt(1), row.getInt(2), row.getDouble(3)))
                    .toList();
        }
    }
}

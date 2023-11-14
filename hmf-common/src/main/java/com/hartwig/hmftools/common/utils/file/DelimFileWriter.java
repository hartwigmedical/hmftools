package com.hartwig.hmftools.common.utils.file;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import org.apache.logging.log4j.util.BiConsumer;


/**
 * Example usage:
 *    new DelimFileWriter().write(fileName, List.of(CHROMOSOME, MEDIAN_RATIO, COUNT), ratios,
 *           (medianRatio, row) -> {
 *               row.set(CHROMOSOME, medianRatio.chromosome());
 *               row.set(MEDIAN_RATIO, medianRatio.medianRatio(), FORMAT);
 *               row.set(COUNT, medianRatio.count()); });
 */
public class DelimFileWriter
{
    // by default, use 4 decimal places for doubles
    private static final NumberFormat sDefaultNumberFormat = new DecimalFormat("#.####", new DecimalFormatSymbols(Locale.ENGLISH));
    private static final String sNullIndicator = "";

    String mDelim = TSV_DELIM;

    public DelimFileWriter()
    {
    }

    public void setDelimiter(String delimiter)
    {
        mDelim = delimiter;
    }

    /**
     * Example usage:
     *    new DelimFileWriter().write(fileName, List.of(CHROMOSOME, MEDIAN_RATIO, COUNT), ratios,
     *           (medianRatio, row) -> {
     *               row.set(CHROMOSOME, medianRatio.chromosome());
     *               row.set(MEDIAN_RATIO, medianRatio.medianRatio(), FORMAT);
     *               row.set(COUNT, medianRatio.count()); });
     *
     * @param filename   path to the file.
     * @param columns    column names of the file.
     * @param objects    iterable of objects to serialised.
     * @param mapper     function to populate a Row given an object of type T.
     */
    public <T> void write(String filename, Iterable<String> columns, Iterable<T> objects, BiConsumer<T, Row> mapper)
    {
        try(BufferedWriter writer = createBufferedWriter(filename))
        {
            write(writer, columns, objects, mapper);
        }
        catch (IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Write to the input buffered writer.
     * @param writer     buffedWriter
     * @param columns    column names of the file.
     * @param objects    iterable of objects to serialised.
     * @param mapper     function to populate a Row given an object of type T.
     */
    public <T> void write(Writer writer, Iterable<String> columns, Iterable<T> objects, BiConsumer<T, Row> mapper)
    {
        try
        {
            // create a map of indices
            Map<String, Integer> columnIndexMap = new HashMap<>();
            int i = 0;
            for (String c : columns)
            {
                if (columnIndexMap.putIfAbsent(c, i++) != null)
                {
                    throw new RuntimeException("duplicate column: " + c);
                }
            }

            writer.write(String.join(mDelim, columns));
            writer.write('\n');

            for (T obj : objects)
            {
                Row row = new Row(columnIndexMap, i); // i is the number of columns
                mapper.accept(obj, row);
                StringJoiner joiner = new StringJoiner(mDelim);
                for(String e : row.mValues)
                {
                    joiner.add(e != null ? e : sNullIndicator);
                }
                writer.write(joiner.toString());
                writer.write('\n');
            }
        }
        catch (IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    // overload that allows using enum as columns
    public <T> void write(String filename, Enum<?>[] columns, Iterable<T> objects, BiConsumer<T, Row> mapper)
    {
        write(filename, Arrays.stream(columns).map(Enum::name).collect(Collectors.toList()), objects, mapper);
    }

    public <T> void write(Writer writer, Enum<?>[] columns, Iterable<T> objects, BiConsumer<T, Row> mapper)
    {
        write(writer, Arrays.stream(columns).map(Enum::name).collect(Collectors.toList()), objects, mapper);
    }

    public static class Row
    {
        private final Map<String, Integer> mColumnIndexMap;
        private final String[] mValues;

        private Row(Map<String, Integer> columnIndexMap, int numColumns)
        {
            mColumnIndexMap = columnIndexMap;
            mValues = new String[numColumns];
        }
        public void set(String column, String value)
        {
            Integer columnIndex = mColumnIndexMap.get(column);
            if (columnIndex == null)
            {
                throw new IllegalArgumentException("invalid column: " + column);
            }
            mValues[columnIndex] = value;
        }
        public void set(String column, int value)
        {
            set(column, Integer.toString(value));
        }

        // we store null as empty string
        public void setNull(String column)
        {
            set(column, sNullIndicator);
        }

        // store bool as 1 and 0
        public void set(String column, boolean value)
        {
            set(column, value ? 1 : 0);
        }

        public void set(String column, char value)
        {
            set(column, String.valueOf(value));
        }

        public void set(String column, byte value)
        {
            set(column, Byte.toString(value));
        }

        public void set(String column, double value)
        {
            set(column, value, sDefaultNumberFormat);
        }

        // row.set("rate", 0.27562, "%.3f");
        public void set(String column, double value, String format)
        {
            set(column, String.format(format, value));
        }

        // row.set("rate", 0.27572, new DecimalFormat("#.####", new DecimalFormatSymbols(Locale.ENGLISH)));
        public void set(String column, double value, NumberFormat format)
        {
            set(column, format.format(value));
        }

        // overloads to allow using enum as column
        public void set(Enum<?> column, String value)
        {
            set(column.name(), value);
        }
        public void set(Enum<?> column, int value)
        {
            set(column.name(), value);
        }
        public void set(Enum<?> column, boolean value)
        {
            set(column.name(), value);
        }
        public void set(Enum<?> column, char value)
        {
            set(column.name(), value);
        }
        public void set(Enum<?> column, byte value)
        {
            set(column.name(), value);
        }
        public void set(Enum<?> column, double value)
        {
            set(column.name(), value);
        }
        public void set(Enum<?> column, double value, String format)
        {
            set(column.name(), value, format);
        }
        public void set(Enum<?> column, double value, NumberFormat format)
        {
            set(column.name(), value, format);
        }
    }
}

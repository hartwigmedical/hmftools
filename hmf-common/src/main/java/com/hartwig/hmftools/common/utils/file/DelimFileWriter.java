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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import org.apache.logging.log4j.util.BiConsumer;

/**
 * Write delimited files such as CSV / TSV.
 *
 * Example usage:
 *
 *    try(DelimFileWriter<Ratio> writer = new DelimFileWriter<>(fileName, List.of(CHROMOSOME, MEDIAN_RATIO, COUNT), ",",
 *           (ratio, row) -> {
 *               row.set(CHROMOSOME, ratio.chromosome());
 *               row.set(MEDIAN_RATIO, ratio.medianRatio(), FORMAT);
 *               row.set(COUNT, ratio.count()); })
 *    {
 *        for(Data d : dataList)
 *        {
 *            writer.writeRow(d);
 *        }
 *    }
 *
 *  There is also a convenience overload to write a file from a collection / iterable of objects:
 *
 *    DelimFileWriter.write(fileName, List.of(CHROMOSOME, MEDIAN_RATIO, COUNT), ratios,
 *           (ratio, row) -> {
 *               row.set(CHROMOSOME, ratio.chromosome());
 *               row.set(MEDIAN_RATIO, ratio.medianRatio(), FORMAT);
 *               row.set(COUNT, ratio.count()); });
 */
public class DelimFileWriter<T> implements AutoCloseable
{
    // by default, use 4 decimal places for doubles
    private static final NumberFormat sDefaultNumberFormat = new DecimalFormat("#.####", new DecimalFormatSymbols(Locale.ENGLISH));
    private static final String sNullIndicator = "null";

    final String mDelimiter;

    final Writer mWriter;

    final List<String> mColumns;

    final BiConsumer<T, Row> mRowEncoder;

    final Map<String, Integer> mColumnIndexMap;

    public DelimFileWriter(String filename, Iterable<String> columns, BiConsumer<T, Row> rowEncoder)
    {
        this(filename, columns, TSV_DELIM, rowEncoder);
    }

    public DelimFileWriter(Writer writer, Iterable<String> columns, BiConsumer<T, Row> rowEncoder)
    {
        this(writer, columns, TSV_DELIM, rowEncoder);
    }

    public DelimFileWriter(String filename, Enum<?>[] columns, BiConsumer<T, Row> rowEncoder)
    {
        this(filename, columns, TSV_DELIM, rowEncoder);
    }

    public DelimFileWriter(Writer writer, Enum<?>[] columns, BiConsumer<T, Row> rowEncoder)
    {
        this(writer, columns, TSV_DELIM, rowEncoder);
    }

    public DelimFileWriter(String filename, Iterable<String> columns, String delimiter, BiConsumer<T, Row> rowEncoder)
    {
        this(createBufferedWriterUnchecked(filename), columns, delimiter, rowEncoder);
    }

    public DelimFileWriter(Writer writer, Iterable<String> columns, String delimiter, BiConsumer<T, Row> rowEncoder)
    {
        mDelimiter = delimiter;
        mWriter = writer;
        mColumns = new ArrayList<>();
        columns.forEach(mColumns::add);
        mRowEncoder = rowEncoder;

        // create a map of indices
        mColumnIndexMap = new HashMap<>();
        int i = 0;
        for(String c : mColumns)
        {
            if(mColumnIndexMap.putIfAbsent(c, i++) != null)
            {
                throw new RuntimeException("duplicate column: " + c);
            }
        }
        try
        {
            mWriter.write(String.join(mDelimiter, mColumns) + '\n');
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    public DelimFileWriter(String filename, Enum<?>[] columns, String delimiter, BiConsumer<T, Row> rowEncoder)
    {
        this(filename, Arrays.stream(columns).map(Enum::name).collect(Collectors.toList()), delimiter, rowEncoder);
    }

    public DelimFileWriter(Writer writer, Enum<?>[] columns, String delimiter, BiConsumer<T, Row> rowEncoder)
    {
        this(writer, Arrays.stream(columns).map(Enum::name).collect(Collectors.toList()), delimiter, rowEncoder);
    }

    // This function is safe to be called from multiple threads.
    // It is user's responsibility to ensure the provided row encoder is thread safe.
    public void writeRow(T obj)
    {
        try
        {
            Row row = new Row(mColumnIndexMap, mColumns.size());
            mRowEncoder.accept(obj, row);
            StringJoiner joiner = new StringJoiner(mDelimiter, "", "\n");
            for(String e : row.mValues)
            {
                joiner.add(e != null ? e : sNullIndicator);
            }

            // NOTE: in order for this function to be thread safe, we must only make a single call to writer
            // otherwise race condition might occur where another thread writes to the file in between
            // writer underneath contains a lock to make sure the write function is thread safe
            mWriter.write(joiner.toString());
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    @Override
    public void close()
    {
        try
        {
            mWriter.close();
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Overload for convenience
     * Example usage:
     *    DelimFileWriter.write(fileName, List.of(CHROMOSOME, MEDIAN_RATIO, COUNT), ratios,
     *           (medianRatio, row) -> {
     *               row.set(CHROMOSOME, medianRatio.chromosome());
     *               row.set(MEDIAN_RATIO, medianRatio.medianRatio(), FORMAT);
     *               row.set(COUNT, medianRatio.count()); });
     *
     * @param filename   path to the file.
     * @param columns    column names of the file.
     * @param objects    iterable of objects to serialised.
     * @param rowEncoder function to populate a Row given an object of type T.
     */
    public static <T> void write(String filename, Iterable<String> columns, Iterable<T> objects, BiConsumer<T, Row> rowEncoder)
    {
        try(BufferedWriter writer = createBufferedWriter(filename))
        {
            write(writer, columns, objects, rowEncoder);
        }
        catch (IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    /**
     * Overload for convenience.
     * Write to the input writer.
     * @param writer     writer
     * @param columns    column names of the file.
     * @param objects    iterable of objects to serialised.
     * @param rowEncoder     function to populate a Row given an object of type T.
     */
    public static <T> void write(Writer writer, Iterable<String> columns, Iterable<T> objects, BiConsumer<T, Row> rowEncoder)
    {
        try(DelimFileWriter<T> delimFileWriter = new DelimFileWriter<>(writer, columns, rowEncoder))
        {
            for(T obj : objects)
            {
                delimFileWriter.writeRow(obj);
            }
        }
    }

    // overload that allows using enum as columns
    public static <T> void write(String filename, Enum<?>[] columns, Iterable<T> objects, BiConsumer<T, Row> rowEncoder)
    {
        write(filename, Arrays.stream(columns).map(Enum::name).collect(Collectors.toList()), objects, rowEncoder);
    }

    public static <T> void write(Writer writer, Enum<?>[] columns, Iterable<T> objects, BiConsumer<T, Row> rowEncoder)
    {
        write(writer, Arrays.stream(columns).map(Enum::name).collect(Collectors.toList()), objects, rowEncoder);
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
        public void set(String key, String value)
        {
            Integer columnIndex = mColumnIndexMap.get(key);
            if (columnIndex == null)
            {
                throw new IllegalArgumentException("invalid column: " + key);
            }
            mValues[columnIndex] = value;
        }
        public void set(String key, int value) { set(key, Integer.toString(value)); }

        // store bool as 1 and 0
        public void set(String key, boolean value) { set(key, value ? 1 : 0); }

        public void set(String key, char value) { set(key, String.valueOf(value)); }

        public void set(String key, byte value) { set(key, Byte.toString(value)); }

        public void set(String key, double value) { set(key, value, sDefaultNumberFormat); }

        //
        //  record.set("rate", 0.27562, "%.3f");
        //
        public void set(String key, double value, String format) { set(key, String.format(format, value)); }

        // record.set("rate", 0.27572, new DecimalFormat("#.####", new DecimalFormatSymbols(Locale.ENGLISH)));
        public void set(String key, double value, NumberFormat format) { set(key, format.format(value)); }

        // overloads to allow using enum as key
        public void set(Enum<?> key, String value) { set(key.name(), value); }
        public void set(Enum<?> key, int value) { set(key.name(), value); }
        public void set(Enum<?> key, boolean value) { set(key.name(), value); }
        public void set(Enum<?> key, char value) { set(key.name(), value); }
        public void set(Enum<?> key, byte value) { set(key.name(), value); }
        public void set(Enum<?> key, double value) { set(key.name(), value); }
        public void set(Enum<?> key, double value, String format) { set(key.name(), value, format); }
        public void set(Enum<?> key, double value, NumberFormat format) { set(key.name(), value, format); }
    }

    // convert to unchecked IO exception to allow usage in streams
    private static BufferedWriter createBufferedWriterUnchecked(String filename)
    {
        try
        {
            return createBufferedWriter(filename);
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }
}

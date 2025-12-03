package com.hartwig.hmftools.common.utils.file;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.Function;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

/**
 * Read a CSV / TSV file. If the file name ends with .gz, it will apply gzip reading.
 * <p>
 * Example usage:
 * try (DelimFileReader reader = new DelimFileReader(filename))
 * {
 * reader.setDelimiter(",");
 * return reader.stream().map(row -> ImmutableMedianRatio.builder()
 * .chromosome(row.get(CHROMOSOME))
 * .medianRatio(row.getDouble(MEDIAN_RATIO))
 * .count(row.getInt(COUNT)).build())
 * .collect(Collectors.toList());
 * }
 * <p>
 * It is also possible to iterate row by row:
 * <p>
 * try (DelimFileReader reader = new DelimFileReader(filename))
 * {
 * for(DelimFileReader.Row row : reader)
 * {
 * String id = row.get("id");
 * String tag = row.getOrNull("tag");
 * int count = row.getInt("count");
 * boolean isValid = row.getBoolean("isValid");
 * Double rate = row.getDoubleOrNull("rate");
 * }
 * }
 * <p>
 * Or use enum instead of String as column identifier:
 * <p>
 * enum Column { id, tag, count, isValid, rate }
 * <p>
 * try (DelimFileReader reader = new DelimFileReader(filename))
 * {
 * for(DelimFileReader.Row row : reader)
 * {
 * String id = row.get(Column.id);
 * String tag = row.getOrNull(Column.tag);
 * int count = row.getInt(Column.count);
 * boolean isValid = row.getBoolean(Column.isValid);
 * Double rate = row.getDoubleOrNull(Column.rate);
 * }
 * }
 */

@SuppressWarnings("unused")
public class DelimFileReader implements Iterable<DelimFileReader.Row>, AutoCloseable
{
    private String mDelim = TSV_DELIM;
    private final BufferedReader mReader;
    private Map<String, Integer> mColumnIndexMap = null;
    private List<String> mColumnNames = null;

    public DelimFileReader(final BufferedReader reader)
    {
        mReader = reader;
    }

    public DelimFileReader(final String filename)
    {
        this(filename, true);
    }

    public DelimFileReader(final String filename, boolean initialiseColumns)
    {
        try
        {
            mReader = createBufferedReader(filename);

            if(initialiseColumns)
            {
                setColumnNames();
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    public void setDelimiter(final String delimiter)
    {
        if(mColumnIndexMap != null)
        {
            throw new IllegalStateException("cannot set delimiter after reading started");
        }
        mDelim = delimiter;
    }

    /**
     * Use the given column names instead of first line of the file.
     *
     * @param columnNames name of columns
     */
    public void setColumnNames(final List<String> columnNames)
    {
        if(mColumnIndexMap != null)
        {
            throw new IllegalStateException("cannot set column names after reading started");
        }
        mColumnIndexMap = new HashMap<>();
        int i = 0;
        for(String columName : columnNames)
        {
            if(mColumnIndexMap.putIfAbsent(columName, i++) != null)
            {
                throw new RuntimeException("duplicate column name: " + columName);
            }
        }
        mColumnNames = Collections.unmodifiableList(columnNames);
    }

    public List<String> getColumnNames()
    {
        if(mColumnIndexMap == null)
        {
            setColumnNames();
        }

        return mColumnNames;
    }

    private void setColumnNames()
    {
        if(mColumnIndexMap == null)
        {
            try
            {
                String line = mReader.readLine();
                setColumnNames(Arrays.asList(line.split(mDelim, -1)));
            }
            catch(IOException e)
            {
                throw new UncheckedIOException(e);
            }
        }
    }

    public boolean hasColumn(final String column)
    {
        return mColumnIndexMap.containsKey(column);
    }

    public boolean hasColumn(final Enum<?> column)
    {
        return hasColumn(column.name());
    }

    // return null if the column is not found
    @Nullable
    public Integer getColumnIndex(final String column)
    {
        return mColumnIndexMap.get(column);
    }

    @Nullable
    public Integer getColumnIndex(final Enum<?> column)
    {
        return getColumnIndex(column.name());
    }

    public Integer getColumnIndex(final Enum<?> column, final String oldNameForColumn)
    {
        Integer columnIndex = getColumnIndex(column.name());
        if(columnIndex == null)
        {
            return getColumnIndex(oldNameForColumn);
        }
        return columnIndex;
    }

    @Override
    public void close()
    {
        try
        {
            mReader.close();
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    @NotNull
    @Override
    public Iterator<Row> iterator()
    {
        // make sure we have got the column names
        setColumnNames();
        return new Iterator<>()
        {
            // almost entirely the same as the one in BufferedReader
            String nextLine = null;

            @Override
            public boolean hasNext()
            {
                if(nextLine != null)
                {
                    return true;
                }
                else
                {
                    try
                    {
                        nextLine = mReader.readLine();
                        return (nextLine != null);
                    }
                    catch(IOException e)
                    {
                        throw new UncheckedIOException(e);
                    }
                }
            }

            @Override
            public Row next()
            {
                if(nextLine != null || hasNext())
                {
                    String line = nextLine;
                    nextLine = null;
                    assert line != null;

                    String[] values = line.split(mDelim, -1);

                    /* revert to in-built split call - performance difference requires more investigation

                    // Split the line on the delimiter. This algorithm is about twice as fast as String.split().
                    String[] values = new String[mColumnIndexMap.size()];
                    int fieldBegin = 0;
                    for(int i = 0; i < values.length; ++i)
                    {
                        int delimIdx = line.indexOf(mDelim, fieldBegin);
                        if(delimIdx < 0)
                        {
                            delimIdx = line.length();
                        }
                        values[i] = line.substring(fieldBegin, delimIdx);
                        fieldBegin = delimIdx + 1;
                    }
                    */

                    return new Row(mColumnIndexMap, values);
                }
                else
                {
                    throw new NoSuchElementException();
                }
            }
        };
    }

    public Stream<Row> stream()
    {
        return StreamSupport.stream(
                Spliterators.spliteratorUnknownSize(
                        iterator(), Spliterator.ORDERED | Spliterator.NONNULL | Spliterator.IMMUTABLE), false);
    }

    public static class Row
    {
        private final Map<String, Integer> mColumnIndexMap;
        private final String[] mValues;

        private Row(final Map<String, Integer> columnIndexMap, final String[] values)
        {
            mColumnIndexMap = columnIndexMap;
            mValues = values;
        }

        public String get(int column) { return getString(column); }
        public String get(String column) { return getString(column); }
        public String get(Enum<?> column) { return getString(column); }
        public String getString(int column) { return getString(column, false); }
        public String getString(String column) { return getString(getColumnIndex(column), false); }
        public String getString(Enum<?> column) { return getString(getColumnIndex(column), false); }
        @Nullable public String getStringOrNull(String column) { return getString(getColumnIndex(column), true); }
        @Nullable public String getStringOrNull(Enum<?> column) { return getString(getColumnIndex(column), true); }

        public byte getByte(String column) { return getByte(getColumnIndex(column), false); }
        public byte getByte(Enum<?> column) { return getByte(getColumnIndex(column), false); }
        @Nullable public Byte getByteOrNull(String column) { return getByte(getColumnIndex(column), true); }
        @Nullable public Byte getByteOrNull(Enum<?> column) { return getByte(getColumnIndex(column), true); }

        public char getChar(String column) { return getChar(getColumnIndex(column), false); }
        public char getChar(Enum<?> column) { return getChar(getColumnIndex(column), false); }
        @Nullable public Character getCharOrNull(String column) { return getChar(getColumnIndex(column), true); }
        @Nullable public Character getCharOrNull(Enum<?> column) { return getChar(getColumnIndex(column), true); }

        public boolean getBoolean(String column) { return getBoolean(getColumnIndex(column), false); }
        public boolean getBoolean(Enum<?> column) { return getBoolean(getColumnIndex(column), false); }
        @Nullable public Boolean getBooleanOrNull(String column) { return getBoolean(getColumnIndex(column), true); }
        @Nullable public Boolean getBooleanOrNull(Enum<?> column) { return getBoolean(getColumnIndex(column), true); }

        public int getInt(int column) { return getInt(column, false); }
        public int getInt(String column) { return getInt(getColumnIndex(column), false); }
        public int getInt(Enum<?> column) { return getInt(getColumnIndex(column), false); }
        @Nullable public Integer getIntOrNull(String column) { return getInt(getColumnIndex(column), true); }
        @Nullable public Integer getIntOrNull(Enum<?> column) { return getInt(getColumnIndex(column), true); }

        public long getLong(String column) { return getLong(getColumnIndex(column), false); }
        public long getLong(Enum<?> column) { return getLong(getColumnIndex(column), false); }
        @Nullable public Long getLongOrNull(String column) { return getLong(getColumnIndex(column), true); }
        @Nullable public Long getLongOrNull(Enum<?> column) { return getLong(getColumnIndex(column), true); }

        public float getFloat(String column) { return getFloat(getColumnIndex(column), false); }
        public float getFloat(Enum<?> column) { return getFloat(getColumnIndex(column), false); }
        @Nullable public Float getFloatOrNull(String column) { return getFloat(getColumnIndex(column), true); }
        @Nullable public Float getFloatOrNull(Enum<?> column) { return getFloat(getColumnIndex(column), true); }

        public double getDouble(int column) { return getDouble(column, false); }
        public double getDouble(String column) { return getDouble(getColumnIndex(column), false); }
        public double getDouble(Enum<?> column) { return getDouble(getColumnIndex(column), false); }
        @Nullable public Double getDoubleOrNull(String column) { return getDouble(getColumnIndex(column), true); }
        @Nullable public Double getDoubleOrNull(Enum<?> column) { return getDouble(getColumnIndex(column), true); }

        private String getString(int column, boolean allowNull)
        {
            return getValue(column, Function.identity(), allowNull);
        }

        private Byte getByte(int column, boolean allowNull)
        {
            return getValue(column, Byte::parseByte, allowNull);
        }

        private Character getChar(int column, boolean allowNull)
        {
            // TODO? should forbid more than 1 char
            return getValue(column, string -> string.charAt(0), allowNull);
        }

        private Boolean getBoolean(int column, boolean allowNull)
        {
            // TODO: should accept only true or false
            return getValue(column, Boolean::parseBoolean, allowNull);
        }

        private Integer getInt(int column, boolean allowNull)
        {
            return getValue(column, Integer::parseInt, allowNull);
        }

        private Long getLong(int column, boolean allowNull)
        {
            return getValue(column, Long::parseLong, allowNull);
        }

        private Float getFloat(int column, boolean allowNull)
        {
            return getValue(column, Float::parseFloat, allowNull);
        }

        private Double getDouble(int column, boolean allowNull)
        {
            return getValue(column, Double::parseDouble, allowNull);
        }

        private <T> T getValue(int column, Function<String, T> parser, boolean allowNull)
        {
            String string = getStringValue(column, allowNull);
            return string == null ? null : parser.apply(string);
        }

        private String getStringValue(int column, boolean allowNull)
        {
            String rawValue = getRawValue(column);
            if (valueIndicatesNull(rawValue))
            {
                if (allowNull)
                {
                    return null;
                }
                else
                {
                    throw new NoSuchElementException(String.format("column %d is null", column));
                }
            }
            else
            {
                return rawValue;
            }
        }

        private int getColumnIndex(Enum<?> column)
        {
            return getColumnIndex(column.name());
        }

        private int getColumnIndex(String column)
        {
            Integer columnIndex = mColumnIndexMap.get(column);
            if(columnIndex == null)
            {
                throw new NoSuchElementException(String.format("column: %s not found", column));
            }
            else
            {
                return columnIndex;
            }
        }

        private String getRawValue(int column)
        {
            return mValues[column];
        }

        private static boolean valueIndicatesNull(String rawValue)
        {
            return rawValue.equals("null") || rawValue.equals("NULL");
        }
    }
}

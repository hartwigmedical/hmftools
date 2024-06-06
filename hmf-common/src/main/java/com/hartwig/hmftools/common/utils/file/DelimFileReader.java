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
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

import org.jetbrains.annotations.Nullable;

/**
 * Read a CSV / TSV file. If the file name ends with .gz, it will apply gzip reading.
 *
 * Example usage:
 *  try (DelimFileReader reader = new DelimFileReader(filename))
 *  {
 *      reader.setDelimiter(",");
 *      return reader.stream().map(row -> ImmutableMedianRatio.builder()
 *                  .chromosome(row.get(CHROMOSOME))
 *                  .medianRatio(row.getDouble(MEDIAN_RATIO))
 *                  .count(row.getInt(COUNT)).build())
 *              .collect(Collectors.toList());
 *  }
 *
 * It is also possible to iterate row by row:
 *
 *  try (DelimFileReader reader = new DelimFileReader(filename))
 *  {
 *      for(DelimFileReader.Row row : reader)
 *      {
 *          String id = row.get("id");
 *          String tag = row.getOrNull("tag");
 *          int count = row.getInt("count");
 *          boolean isValid = row.getBoolean("isValid");
 *          Double rate = row.getDoubleOrNull("rate");
 *      }
 *  }
 *
 *  Or use enum instead of String as column identifier:
 *
 *  enum Column { id, tag, count, isValid, rate }
 *
 *  try (DelimFileReader reader = new DelimFileReader(filename))
 *  {
 *      for(DelimFileReader.Row row : reader)
 *      {
 *          String id = row.get(Column.id);
 *          String tag = row.getOrNull(Column.tag);
 *          int count = row.getInt(Column.count);
 *          boolean isValid = row.getBoolean(Column.isValid);
 *          Double rate = row.getDoubleOrNull(Column.rate);
 *      }
 *  }
 */
public class DelimFileReader implements Iterable<DelimFileReader.Row>, AutoCloseable
{
    private String mDelim = TSV_DELIM;
    private final BufferedReader mReader;
    private Map<String,Integer> mColumnIndexMap = null;
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
            setColumnNames();

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

    public boolean hasColumn(final String column) { return mColumnIndexMap.containsKey(column); }
    public boolean hasColumn(final Enum<?> column) { return hasColumn(column.name()); }

    // return null if the column is not found
    @Nullable
    public Integer getColumnIndex(final String column) { return mColumnIndexMap.get(column); }
    @Nullable
    public Integer getColumnIndex(final Enum<?> column) { return getColumnIndex(column.name()); }

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

    @Override
    public Iterator<Row> iterator()
    {
        // make sure we have got the column names
        getColumnNames();
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

        public boolean isNull(final String column)
        {
            return valueIndicatesNull(parseRawValue(column));
        }

        public String get(final String column)
        {
            String v = parseRawValue(column);
            if(valueIndicatesNull(v))
            {
                throw new NoSuchElementException();
            }
            return v;
        }

        public @Nullable String getOrNull(final String column)
        {
            String v = parseRawValue(column);
            if(valueIndicatesNull(v))
                return null;
            return v;
        }

        public int getInt(final String column)
        {
            return Integer.parseInt(get(column));
        }

        public @Nullable Integer getIntOrNull(final String column)
        {
            String v = getOrNull(column);
            return v == null ? null : Integer.parseInt(v);
        }

        // store boolean as 1 and 0
        public boolean getBoolean(final String column)
        {
            return getInt(column) != 0;
        }

        public @Nullable Boolean getBooleanOrNull(final String column)
        {
            Integer v = getIntOrNull(column);
            return v == null ? null : v != 0;
        }

        public char getChar(final String column)
        {
            return get(column).charAt(0);
        }

        public @Nullable Character getCharOrNull(final String column)
        {
            String v = getOrNull(column);
            return v == null ? null : v.charAt(0);
        }

        public byte getByte(final String column)
        {
            return Byte.parseByte(get(column));
        }

        public @Nullable Byte getByteOrNull(final String column)
        {
            String v = getOrNull(column);
            return v == null ? null : Byte.parseByte(v);
        }

        public double getDouble(final String column)
        {
            return Double.parseDouble(get(column));
        }

        public @Nullable Double getDoubleOrNull(final String column)
        {
            String v = getOrNull(column);
            return v == null ? null : Double.parseDouble(v);
        }

        public long getLong(final String column)
        {
            return Long.parseLong(get(column));
        }

        public @Nullable Long getLongOrNull(final String column)
        {
            String v = getOrNull(column);
            return v == null ? null : Long.parseLong(v);
        }

        // overloads that allow using enum as column

        public boolean isNull(final Enum<?> column)
        {
            return isNull(column.name());
        }

        public String get(final Enum<?> column)
        {
            return get(column.name());
        }
        public @Nullable String getOrNull(final Enum<?> column)
        {
            return getOrNull(column.name());
        }

        public int getInt(final Enum<?> column)
        {
            return getInt(column.name());
        }
        public @Nullable Integer getIntOrNull(final Enum<?> column)
        {
            return getIntOrNull(column.name());
        }

        public boolean getBoolean(final Enum<?> column)
        {
            return getBoolean(column.name());
        }
        public @Nullable Boolean getBooleanOrNull(final Enum<?> column)
        {
            return getBooleanOrNull(column.name());
        }

        public char getChar(final Enum<?> column)
        {
            return getChar(column.name());
        }
        public @Nullable Character getCharOrNull(final Enum<?> column)
        {
            return getCharOrNull(column.name());
        }

        public byte getByte(final Enum<?> column)
        {
            return getByte(column.name());
        }
        public @Nullable Byte getByteOrNull(final Enum<?> column)
        {
            return getByteOrNull(column.name());
        }

        public double getDouble(final Enum<?> column)
        {
            return getDouble(column.name());
        }
        public @Nullable Double getDoubleOrNull(final Enum<?> column)
        {
            return getDoubleOrNull(column.name());
        }

        public long getLong(final Enum<?> column)
        {
            return getLong(column.name());
        }
        public @Nullable Long getLongOrNull(final Enum<?> column)
        {
            return getLongOrNull(column.name());
        }

        // overloads that get by column index, provided for speed reasons.
        // Onus is on caller to make sure the indices are valid
        // Or IndexOutOfBoundsException will be thrown
        public String get(final int columnIndex)
        {
            String v = mValues[columnIndex];
            if(valueIndicatesNull(v))
            {
                throw new NoSuchElementException();
            }
            return v;
        }

        public @Nullable String getOrNull(final int columnIndex)
        {
            String v = mValues[columnIndex];
            if(valueIndicatesNull(v))
                return null;
            return v;
        }

        public int getInt(final int columnIndex)
        {
            return Integer.parseInt(get(columnIndex));
        }

        public @Nullable Integer getIntOrNull(final int columnIndex)
        {
            String v = getOrNull(columnIndex);
            return v == null ? null : Integer.parseInt(v);
        }

        // store boolean as 1 and 0
        public boolean getBoolean(final int columnIndex)
        {
            return getInt(columnIndex) != 0;
        }

        public @Nullable Boolean getBooleanOrNull(final int columnIndex)
        {
            Integer v = getIntOrNull(columnIndex);
            return v == null ? null : v != 0;
        }

        public char getChar(final int columnIndex)
        {
            return get(columnIndex).charAt(0);
        }

        public @Nullable Character getCharOrNull(final int columnIndex)
        {
            String v = getOrNull(columnIndex);
            return v == null ? null : v.charAt(0);
        }

        public byte getByte(final int columnIndex)
        {
            return Byte.parseByte(get(columnIndex));
        }

        public @Nullable Byte getByteOrNull(final int columnIndex)
        {
            String v = getOrNull(columnIndex);
            return v == null ? null : Byte.parseByte(v);
        }

        public double getDouble(final int columnIndex)
        {
            return Double.parseDouble(get(columnIndex));
        }

        public @Nullable Double getDoubleOrNull(final int columnIndex)
        {
            String v = getOrNull(columnIndex);
            return v == null ? null : Double.parseDouble(v);
        }

        public long getLong(final int columnIndex)
        {
            return Long.parseLong(get(columnIndex));
        }

        public @Nullable Long getLongOrNull(final int columnIndex)
        {
            String v = getOrNull(columnIndex);
            return v == null ? null : Long.parseLong(v);
        }

        // get the value that is stored in the row
        private String parseRawValue(String column)
        {
            Integer index = mColumnIndexMap.get(column);
            if(index == null)
            {
                throw new NoSuchElementException(String.format("column: %s not found", column));
            }
            return mValues[index];
        }

        private static boolean valueIndicatesNull(String rawValue)
        {
            return rawValue.equals("null") || rawValue.equals("NULL");
        }
    }
}

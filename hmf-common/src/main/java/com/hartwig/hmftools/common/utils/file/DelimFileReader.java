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

import org.jetbrains.annotations.NotNull;
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
 * It is also possible to iterate row by row
 *
 *  try (DelimFileReader reader = new DelimFileReader(filename))
 *  {
 *      for(DelimFileReader.Row row : reader)
 *      {
 *          String id = row.get("id");
 *          int count = row.get("count");
 *          boolean isValid = row.getBoolean("isValid");
 *          Double rate = row.getDoubleOrNull("rate");
 *      }
 *  }
 */
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
        try
        {
            mReader = createBufferedReader(filename);
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
        return mColumnNames;
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

        private Row(Map<String, Integer> columnIndexMap, String[] values)
        {
            mColumnIndexMap = columnIndexMap;
            mValues = values;
        }

        public boolean isNull(String column)
        {
            return valueIndicatesNull(parseRawValue(column));
        }

        public String get(String column)
        {
            String v = parseRawValue(column);
            if(valueIndicatesNull(v))
            {
                throw new NoSuchElementException();
            }
            return v;
        }

        public @Nullable String getOrNull(String column)
        {
            String v = parseRawValue(column);
            if(valueIndicatesNull(v))
                return null;
            return v;
        }

        public int getInt(String column)
        {
            return Integer.parseInt(get(column));
        }

        public @Nullable Integer getIntOrNull(String column)
        {
            String v = getOrNull(column);
            return v == null ? null : Integer.parseInt(v);
        }

        // store boolean as 1 and 0
        public boolean getBoolean(String column)
        {
            return getInt(column) != 0;
        }

        public @Nullable Boolean getBooleanOrNull(String column)
        {
            Integer v = getIntOrNull(column);
            return v == null ? null : v != 0;
        }

        public char getChar(String column)
        {
            return get(column).charAt(0);
        }

        public @Nullable Character getCharOrNull(String column)
        {
            String v = getOrNull(column);
            return v == null ? null : v.charAt(0);
        }

        public byte getByte(String column)
        {
            return Byte.parseByte(get(column));
        }

        public @Nullable Byte getByteOrNull(String column)
        {
            String v = getOrNull(column);
            return v == null ? null : Byte.parseByte(v);
        }

        public double getDouble(String column)
        {
            return Double.parseDouble(get(column));
        }

        public @Nullable Double getDoubleOrNull(String column)
        {
            String v = getOrNull(column);
            return v == null ? null : Double.parseDouble(v);
        }

        public long getLong(String column)
        {
            return Long.parseLong(get(column));
        }

        public @Nullable Long getLongOrNull(String column)
        {
            String v = getOrNull(column);
            return v == null ? null : Long.parseLong(v);
        }

        // overloads that allow using enum as column
        public boolean isNull(Enum<?> column)
        {
            return isNull(column.name());
        }

        public String get(Enum<?> column)
        {
            return get(column.name());
        }
        public @Nullable String getOrNull(Enum<?> column)
        {
            return getOrNull(column.name());
        }

        public int getInt(Enum<?> column)
        {
            return getInt(column.name());
        }
        public @Nullable Integer getIntOrNull(Enum<?> column)
        {
            return getIntOrNull(column.name());
        }

        public boolean getBoolean(Enum<?> column)
        {
            return getBoolean(column.name());
        }
        public @Nullable Boolean getBooleanOrNull(Enum<?> column)
        {
            return getBooleanOrNull(column.name());
        }

        public char getChar(Enum<?> column)
        {
            return getChar(column.name());
        }
        public @Nullable Character getCharOrNull(Enum<?> column)
        {
            return getCharOrNull(column.name());
        }

        public byte getByte(Enum<?> column)
        {
            return getByte(column.name());
        }
        public @Nullable Byte getByteOrNull(Enum<?> column)
        {
            return getByteOrNull(column.name());
        }

        public double getDouble(Enum<?> column)
        {
            return getDouble(column.name());
        }
        public @Nullable Double getDoubleOrNull(Enum<?> column)
        {
            return getDoubleOrNull(column.name());
        }

        public long getLong(Enum<?> column)
        {
            return getLong(column.name());
        }
        public @Nullable Long getLongOrNull(Enum<?> column)
        {
            return getLongOrNull(column.name());
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
            return rawValue.isEmpty() || rawValue.equals("null") || rawValue.equals("NULL");
        }
    }
}

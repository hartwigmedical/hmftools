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
 *     try (DelimFileReader reader = new DelimFileReader(filename))
 *     {
 *         reader.setDelimiter(",");
 *         return reader.stream().map(row -> ImmutableMedianRatio.builder()
 *                     .chromosome(row.get(CHROMOSOME))
 *                     .medianRatio(row.getDouble(MEDIAN_RATIO))
 *                     .count(row.getInt(COUNT)).build())
 *                 .collect(Collectors.toList());
 *     }
 *
 * It is also possible to iterate row by row
 *
 *  try (DelimFileReader reader = new DelimFileReader(filename))
 *  {
 *      for(DelimFileReader.Row row : reader)
 *      {
 *          int id = row.getInt("id");
 *      }
 *  }
 */
public class DelimFileReader implements Iterable<DelimFileReader.Row>, AutoCloseable
{
    private String mDelim = TSV_DELIM;
    private final BufferedReader mReader;
    private Map<String, Integer> mColumnIndexMap = null;
    private List<String> mColumnNames = null;

    public DelimFileReader(@NotNull BufferedReader reader)
    {
        mReader = reader;
    }

    public DelimFileReader(@NotNull String filename)
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

    public void setDelimiter(@NotNull String delimiter)
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
    public void setColumnNames(@NotNull List<String> columnNames)
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

    @NotNull
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

    @NotNull
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

    @NotNull
    public Stream<Row> stream()
    {
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(iterator(), Spliterator.ORDERED | Spliterator.NONNULL | Spliterator.IMMUTABLE), false);
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

        public boolean isNull(String key)
        {
            return valueIndicatesNull(parseRawValue(key));
        }

        public @Nullable String get(String key)
        {
            String v = parseRawValue(key);
            if(valueIndicatesNull(v))
                return null;
            return v;
        }

        public @Nullable Integer getInt(String key)
        {
            String v = get(key);
            if(v == null)
                return null;
            return Integer.parseInt(v);
        }

        // store boolean as 1 and 0
        public @Nullable Boolean getBoolean(String key)
        {
            Integer v = getInt(key);
            if(v == null)
                return null;
            return v == 1;
        }

        public @Nullable Character getChar(String key)
        {
            String v = get(key);
            if(v == null)
                return null;
            return v.charAt(0);
        }

        public @Nullable Byte getByte(String key)
        {
            String v = get(key);
            if(v == null)
                return null;
            return Byte.parseByte(v);
        }

        public @Nullable Double getDouble(String key)
        {
            String v = get(key);
            if(v == null)
                return null;
            return Double.parseDouble(v);
        }

        // overloads that allow using enum as key
        public boolean isNull(Enum<?> key)
        {
            return isNull(key.name());
        }

        public @Nullable String get(Enum<?> key)
        {
            return get(key.name());
        }

        public @Nullable Integer getInt(Enum<?> key)
        {
            return getInt(key.name());
        }

        public @Nullable Boolean getBoolean(Enum<?> key)
        {
            return getBoolean(key.name());
        }

        public @Nullable Character getChar(Enum<?> key)
        {
            return getChar(key.name());
        }

        public @Nullable Byte getByte(Enum<?> key)
        {
            return getByte(key.name());
        }

        public @Nullable Double getDouble(Enum<?> key)
        {
            return getDouble(key.name());
        }

        // get the value that is stored in the file
        private @NotNull String parseRawValue(String key)
        {
            Integer index = mColumnIndexMap.get(key);
            if(index == null)
            {
                throw new RuntimeException(String.format("column: %s not found", key));
            }
            return mValues[index];
        }

        private static boolean valueIndicatesNull(@NotNull String rawValue)
        {
            return rawValue.isEmpty() || rawValue.equals("null") || rawValue.equals("NULL");
        }
    }
}

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
        catch (IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    public void setDelimiter(@NotNull String delimiter)
    {
        if (mColumnIndexMap != null)
        {
            throw new IllegalStateException("cannot set delimiter after reading started");
        }
        mDelim = delimiter;
    }

    /**
     * Use the given column names instead of first line of the file.
     * @param columnNames  name of columns
     */
    public void setColumnNames(@NotNull List<String> columnNames)
    {
        if (mColumnIndexMap != null)
        {
            throw new IllegalStateException("cannot set column names after reading started");
        }
        mColumnIndexMap = new HashMap<>();
        int i = 0;
        for (String columName : columnNames)
        {
            if (mColumnIndexMap.putIfAbsent(columName, i++) != null)
            {
                throw new RuntimeException("duplicate column name: " + columName);
            }
        }
        mColumnNames = Collections.unmodifiableList(columnNames);
    }

    @NotNull
    public List<String> getColumnNames()
    {
        if (mColumnIndexMap == null)
        {
            try
            {
                String line = mReader.readLine();
                setColumnNames(Arrays.asList(line.split(mDelim, -1)));
            }
            catch (IOException e)
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
        catch (IOException e)
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
                if (nextLine != null)
                {
                    return true;
                }
                else
                {
                    try {
                        nextLine = mReader.readLine();
                        return (nextLine != null);
                    }
                    catch (IOException e) {
                        throw new UncheckedIOException(e);
                    }
                }
            }

            @Override
            public Row next()
            {
                if (nextLine != null || hasNext())
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
        return StreamSupport.stream(Spliterators.spliteratorUnknownSize(
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
        public boolean isNull(String key)
        {
            String v = get(key);
            return !v.isEmpty() && !v.equals("null") && !v.equals("NULL");
        }
        public String get(String key)
        {
            Integer index = mColumnIndexMap.get(key);
            if (index == null)
            {
                throw new RuntimeException(String.format("column: %s not found", key));
            }
            return mValues[index];
        }
        public Integer getInt(String key)
        {
            // TODO: handle null
            return Integer.parseInt(get(key));
        }
        public Double getDouble(String key)
        {
            // TODO: handle null
            return Double.parseDouble(get(key));
        }
        public Boolean getBool(String key)
        {
            // TODO: handle null
            return Boolean.parseBoolean(get(key));
        }

        // overloads that allow using enum as key
        public boolean isNull(Enum<?> key)
        {
            return isNull(key.name());
        }
        public String get(Enum<?> key)
        {
            return get(key.name());
        }
        public Integer getInt(Enum<?> key)
        {
            return getInt(key.name());
        }
        public Double getDouble(Enum<?> key)
        {
            return getDouble(key.name());
        }
        public Boolean getBool(Enum<?> key)
        {
            return getBool(key.name());
        }
    }
}

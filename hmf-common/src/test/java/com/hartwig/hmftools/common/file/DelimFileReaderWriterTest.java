package com.hartwig.hmftools.common.file;

import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.junit.Test;

public class DelimFileReaderWriterTest
{
    static class Data
    {
        public int count;
        public String name;
        public double rate;

        public Data(int count, String name, double rate)
        {
            this.count = count;
            this.name = name;
            this.rate = rate;
        }

        public void assertValuesEqual(int count, String name, double rate)
        {
            assertEquals(this.count, count);
            assertEquals(this.name, name);
            assertEquals(this.rate, rate, 1e-5);
        }
    }

    @Test
    public void testWriter()
    {
        StringWriter stringWriter = new StringWriter();

        List<Data> dataList = new ArrayList<>();
        dataList.add(new Data(10, "Susan", 0.7));
        dataList.add(new Data(29, "John", 0.8));

        try (BufferedWriter bufferedWriter = new BufferedWriter(stringWriter))
        {
            DelimFileWriter<Data> writer = new DelimFileWriter<>(bufferedWriter, List.of("count", "name", "rate"), ",", ((data, row) ->
                {
                    row.set("name", data.name);
                    row.set("count", data.count);
                    row.set("rate", data.rate);
                }));
            dataList.forEach(writer::writeRow);
        }
        catch (IOException e)
        {
            throw new RuntimeException(e);
        }

        assertEquals("count,name,rate\n" + "10,Susan,0.7\n" + "29,John,0.8\n", stringWriter.toString());
    }

    @Test
    public void testReader()
    {
        String csvContent = "count,name,rate\n" +
                            "10,Susan,0.7\n" +
                            "29,John,0.8\n";

        try(DelimFileReader reader = new DelimFileReader(new BufferedReader(new StringReader(csvContent))))
        {
            reader.setDelimiter(",");
            List<Data> dataList = reader.stream()
                    .map(row -> new Data(row.getInt("count"), row.get("name"), row.getDouble("rate")))
                    .collect(Collectors.toList());

            assertEquals(2, dataList.size());
            dataList.get(0).assertValuesEqual(10, "Susan", 0.7);
            dataList.get(1).assertValuesEqual(29, "John", 0.8);
        }

        // test reading without stream
        try(DelimFileReader reader = new DelimFileReader(new BufferedReader(new StringReader(csvContent))))
        {
            reader.setDelimiter(",");
            List<Data> dataList = new ArrayList<>();
            for(DelimFileReader.Row row : reader)
            {
                dataList.add(new Data(row.getInt(0), row.get(1), row.getDouble(2)));
            }

            assertEquals(2, dataList.size());
            dataList.get(0).assertValuesEqual(10, "Susan", 0.7);
            dataList.get(1).assertValuesEqual(29, "John", 0.8);
        }
    }
}

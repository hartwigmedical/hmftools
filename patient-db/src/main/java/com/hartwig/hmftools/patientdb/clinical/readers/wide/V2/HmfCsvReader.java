package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

public final class HmfCsvReader
{
    private HmfCsvReader()
    {
    }

    public static List<CsvEntry> read(String pathToCsv) throws IOException
    {
        InputStream file = Thread.currentThread().getContextClassLoader().getResourceAsStream(pathToCsv);
        if(file == null)
        {
            throw new FileNotFoundException(String.format("File '%s' not found.", pathToCsv));
        }
        return read(file);
    }

    public static List<CsvEntry> read(InputStream inputStream) throws IOException
    {
        BufferedReader csvReader = new BufferedReader(new InputStreamReader(inputStream));
        String headerLine = csvReader.readLine();
        if(headerLine == null || headerLine.isBlank())
        {
            return List.of();
        }
        String[] headers = headerLine.split(",");
        return csvReader.lines().filter(line -> !line.isBlank())
                .map(HmfCsvReader::removeQuoteCharacters)
                .map(line -> line.split(","))
                .map(columns -> new CsvEntry(headers, columns))
                .collect(Collectors.toList());
    }

    private static String removeQuoteCharacters(String s)
    {
        return s.replace("\"", "");
    }

}

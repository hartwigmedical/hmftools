package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

import org.jooq.tools.csv.CSVReader;

public final class HmfCsvReader
{
    private HmfCsvReader()
    {
    }

    public static List<CsvEntry> read(String pathToCsv) throws IOException
    {
        InputStream file = Thread.currentThread().getContextClassLoader().getResourceAsStream(pathToCsv);
        if (file == null) {
            throw new FileNotFoundException(String.format("File '%s' not found.", pathToCsv));
        }
        return read(file);
    }

    public static List<CsvEntry> read(InputStream inputStream) throws IOException
    {
        CSVReader csvReader = new CSVReader(new InputStreamReader(inputStream));
        String[] headers = csvReader.readNext();
        return csvReader.readAll().stream().map(values -> new CsvEntry(headers, values)).collect(Collectors.toList());
    }

}

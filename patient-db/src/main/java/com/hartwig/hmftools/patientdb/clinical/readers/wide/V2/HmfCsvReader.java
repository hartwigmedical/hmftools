package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.jooq.tools.csv.CSVReader;

public final class HmfCsvReader
{
    private HmfCsvReader()
    {
    }

    public static List<CsvEntry> read(String pathToCsv) throws IOException
    {
        List<CsvEntry> result = new ArrayList<>();
        CSVReader csvReader = new CSVReader(new FileReader(pathToCsv));
        String[] headers = csvReader.readNext();

        while(csvReader.hasNext())
        {
            String[] values = csvReader.readNext();
            CsvEntry entry = new CsvEntry(headers, values);
            result.add(entry);
        }
        return result;
    }

}

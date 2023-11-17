package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CsvEntry
{
    private final Map<String, String> entry;

    public CsvEntry(String[] headers, String[] values)
    {
        if(headers.length != values.length)
        {
            throw new IllegalArgumentException(String.format("Amount of headers was '%d' and amount of values was '%d' but they should be equal.", headers.length, values.length));
        }
        entry = IntStream.range(0, headers.length).boxed().collect(Collectors.toMap(i -> headers[i], i -> values[i]));
    }

    public Optional<String> get(String header)
    {
        if(!entry.containsKey(header))
        {
            throw new NoSuchElementException("Header with name '%s' does not exist.");
        }
        var result = entry.get(header);
        if(result.equals("NA"))
        {
            return Optional.empty();
        }
        return Optional.of(result);
    }
}

package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;

import org.junit.Test;

public class HmfCsvReaderTest
{
    private static final String DUMMY_CSV_PATH = "readers/dummy_csv_file.csv";

    @Test
    public void canReadDummyFile() throws IOException
    {
        List<CsvEntry> entries = HmfCsvReader.read(DUMMY_CSV_PATH);

        assertEquals(2, entries.size());

        var firstEntry = entries.get(0);
        var firstEntryValue0 = firstEntry.get("header1").orElseThrow();
        var firstEntryValue1 = firstEntry.get("header2").orElseThrow();
        var firstEntryValue2 = firstEntry.get("header3").orElseThrow();

        assertEquals("0", firstEntryValue0);
        assertEquals("foo", firstEntryValue1);
        assertEquals("bar", firstEntryValue2);

        var secondEntry = entries.get(1);
        var secondEntryValue0 = secondEntry.get("header1").orElseThrow();
        var secondEntryValue1 = secondEntry.get("header2").orElseThrow();
        var secondEntryValue2 = secondEntry.get("header3").orElseThrow();

        assertEquals("1", secondEntryValue0);
        assertEquals("baz", secondEntryValue1);
        assertEquals("qux", secondEntryValue2);
    }

    @Test
    public void nonExistentFileThrows() throws IOException
    {
        assertThrows(FileNotFoundException.class, () -> HmfCsvReader.read("NonExistentCsvPath.csv"));
    }

    /**
     * A bug occurred where an exception would be thrown if a csv file was parsed that had a newline at the end.
     * This test should catch that bug.
     */
    @Test
    public void testHeaderOnlyWithNewlineDoesNotThrow() throws IOException
    {
        List<CsvEntry> result = HmfCsvReader.read("readers/file_with_newlines_at_end.csv");
        assertEquals(1, result.size());
    }

    @Test
    public void testCanReadFileWithQuotedFields() throws IOException
    {
        List<CsvEntry> result = HmfCsvReader.read("readers/file_with_quoted_fields.csv");
        assertEquals(1, result.size());

        var entry = result.get(0);
        var column = entry.get("header1").orElseThrow();
        assertEquals("value", column);
    }
}

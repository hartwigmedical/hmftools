package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.NoSuchElementException;

import org.junit.Test;

public class CsvEntryTest
{

    @Test
    public void testInitialisationEmptyDoesNotThrow()
    {
        try
        {
            new CsvEntry(new String[] {}, new String[] {});
        }
        catch(Exception e)
        {
            fail("Unexpected exception was thrown: " + e);
        }
    }

    @Test
    public void testInitWithDifferentLengthShouldThrow()
    {
        var headers = new String[] { "header1" };
        var valuesWithDifferentLength = new String[] { "value1", "extra-value" };

        assertThrows(IllegalArgumentException.class, () -> new CsvEntry(headers, valuesWithDifferentLength));
    }

    @Test
    public void testSimpleEntry()
    {
        var headers = new String[] { "header1", "header2" };
        var values = new String[] { "value1", "value2" };
        var entry = new CsvEntry(headers, values);

        var resolvedValue = entry.get("header1").orElseThrow();
        assertEquals("value1", resolvedValue);
    }

    @Test
    public void testNaReturnsEmptyOptional()
    {
        var headers = new String[] { "header1" };
        var values = new String[] { "NA" };
        var entry = new CsvEntry(headers, values);

        assertTrue(entry.get("header1").isEmpty());
    }

    @Test
    public void testAccessingNonExistentHeaderThrows()
    {
        var headers = new String[] { "header1" };
        var values = new String[] { "value1" };
        var entry = new CsvEntry(headers, values);

        assertThrows(NoSuchElementException.class, () -> entry.get("non existent header"));
    }

}

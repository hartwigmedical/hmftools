package com.hartwig.hmftools.serve.sources.hartwig;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class HartwigFileReaderTest {

    private static final String EXAMPLE_TSV = Resources.getResource("hartwig/example.tsv").getPath();

    @Test
    public void canReadHartwigFile() throws IOException {
        List<HartwigEntry> entries = HartwigFileReader.read(EXAMPLE_TSV);

        assertEquals(3, entries.size());

        assertEquals("5", entries.get(0).chromosome());
        assertEquals(1295228, entries.get(0).position());
        assertEquals("G", entries.get(0).ref());
        assertEquals("A", entries.get(0).alt());
        assertEquals("TERT", entries.get(0).gene());
        assertEquals("ENST00000310581", entries.get(0).transcript());
        assertEquals(Strings.EMPTY, entries.get(0).proteinAnnotation());

        assertEquals("7", entries.get(1).chromosome());
        assertEquals(140453136, entries.get(1).position());
        assertEquals("A", entries.get(1).ref());
        assertEquals("T", entries.get(1).alt());
        assertEquals("BRAF", entries.get(1).gene());
        assertEquals("ENST00000288602", entries.get(1).transcript());
        assertEquals("V600E", entries.get(1).proteinAnnotation());

        assertEquals("7", entries.get(2).chromosome());
        assertEquals(140453136, entries.get(2).position());
        assertEquals("A", entries.get(2).ref());
        assertEquals("T", entries.get(2).alt());
        assertEquals("BRAF", entries.get(2).gene());
        assertEquals("ENST00000288602", entries.get(2).transcript());
        assertEquals(Strings.EMPTY, entries.get(2).proteinAnnotation());
    }
}
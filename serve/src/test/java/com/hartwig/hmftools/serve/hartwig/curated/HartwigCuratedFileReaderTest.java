package com.hartwig.hmftools.serve.hartwig.curated;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class HartwigCuratedFileReaderTest {

    private static final String EXAMPLE_TSV = Resources.getResource("hartwig/curated/example.tsv").getPath();

    @Test
    public void canReadHartwigCohortFile() throws IOException {
        List<HartwigCuratedEntry> entries = HartwigCuratedFileReader.readCuratedFile(EXAMPLE_TSV);

        assertEquals(2, entries.size());

        assertEquals("1", entries.get(0).chromosome());
        assertEquals(226252135, entries.get(0).position());
        assertEquals("A", entries.get(0).ref());
        assertEquals("T", entries.get(0).alt());
        assertEquals("H3F3A", entries.get(0).gene());
        assertEquals("ENST00000366813", entries.get(0).transcript());
        assertEquals("K28M", entries.get(0).proteinAnnotation());

        assertEquals("1", entries.get(1).chromosome());
        assertEquals(149812647, entries.get(1).position());
        assertEquals("T", entries.get(1).ref());
        assertEquals("A", entries.get(1).alt());
        assertEquals("HIST2H3C", entries.get(1).gene());
        assertEquals("ENST00000369158", entries.get(1).transcript());
        assertEquals("K28M", entries.get(1).proteinAnnotation());
    }

}
package com.hartwig.hmftools.fastqstats;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class FastqTrackerTest {

    private static final String BASE_PATH = Resources.getResource("fastq").getPath();

    @Test
    public void computesCorrectPercentageFor100PercentFile() throws IOException {
        final String path = BASE_PATH + File.separator + "q30-100.fastq";
        final FileInputStream fin = new FileInputStream(path);
        final FastqReader fr = new FastqReader(fin);
        final FastqData data = fr.read();
        assertEquals(100, data.yield());
        assertEquals(100.0, data.q30Percentage(), 0.00001);
    }

    @Test
    public void computesCorrectPercentageFor0PercentFile() throws IOException {
        final String path = BASE_PATH + File.separator + "q30-0.fastq";
        final FileInputStream fin = new FileInputStream(path);
        final FastqReader fr = new FastqReader(fin);
        final FastqData data = fr.read();
        assertEquals(100, data.yield());
        assertEquals(0.0, data.q30Percentage(), 0.00001);
    }

    @Test
    public void computesCorrectPercentageFor10PercentFile() throws IOException {
        final String path = BASE_PATH + File.separator + "q30-10.fastq";
        final FileInputStream fin = new FileInputStream(path);
        final FastqReader fr = new FastqReader(fin);
        final FastqData data = fr.read();
        assertEquals(100, data.yield());
        assertEquals(10.0, data.q30Percentage(), 0.00001);
    }
}

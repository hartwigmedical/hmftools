package com.hartwig.hmftools.fastqstats;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.net.URL;

import com.google.common.io.Resources;

import org.junit.Test;

public class FastqTrackerTest {
    @Test
    public void computesCorrectPercentageFor100PercentFile() throws IOException {
        final URL path = Resources.getResource("fastq/q30-100.fastq");
        final FileInputStream fin = new FileInputStream(new File(path.getPath()));
        final FastqReader fr = new FastqReader(fin);
        final FastqData data = fr.read();
        final double q30Percentage = data.getQ30() * 100.0 / data.getYield();
        assertEquals(100, data.getYield());
        assertEquals(100.0, q30Percentage, 0.00001);
    }

    @Test
    public void computesCorrectPercentageFor0PercentFile() throws IOException {
        final URL path = Resources.getResource("fastq/q30-0.fastq");
        final FileInputStream fin = new FileInputStream(new File(path.getPath()));
        final FastqReader fr = new FastqReader(fin);
        final FastqData data = fr.read();
        final double q30Percentage = data.getQ30() * 100.0 / data.getYield();
        assertEquals(100, data.getYield());
        assertEquals(0.0, q30Percentage, 0.00001);
    }

    @Test
    public void computesCorrectPercentageFor10PercentFile() throws IOException {
        final URL path = Resources.getResource("fastq/q30-10.fastq");
        final FileInputStream fin = new FileInputStream(new File(path.getPath()));
        final FastqReader fr = new FastqReader(fin);
        final FastqData data = fr.read();
        final double q30Percentage = data.getQ30() * 100.0 / data.getYield();
        assertEquals(100, data.getYield());
        assertEquals(10.0, q30Percentage, 0.00001);
    }
}

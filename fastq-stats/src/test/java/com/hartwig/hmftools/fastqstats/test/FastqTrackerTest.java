package com.hartwig.hmftools.fastqstats.test;

import org.junit.Test;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.fastqstats.FastqData;
import com.hartwig.hmftools.fastqstats.FastqReader;

public class FastqTrackerTest {
    @Test
    public void computesCorrectPercentageFor100PercentFile() throws IOException{
        FileInputStream fin = new FileInputStream(new File("./src/test/resources/fastq/q30-100.fastq"));
        FastqReader fr = new FastqReader(fin);
        FastqData data = fr.read();
        double q30Percentage = data.getQ30() * 100.0 / data.getYield();
        assertEquals(100, data.getYield());
        assertEquals(100.0, q30Percentage, 0.00001);
    }

    @Test
    public void computesCorrectPercentageFor0PercentFile() throws IOException{
        FileInputStream fin = new FileInputStream(new File("./src/test/resources/fastq/q30-0.fastq"));
        FastqReader fr = new FastqReader(fin);
        FastqData data = fr.read();
        double q30Percentage = data.getQ30() * 100.0 / data.getYield();
        assertEquals(100, data.getYield());
        assertEquals(0.0, q30Percentage,0.00001);
    }

    @Test
    public void computesCorrectPercentageFor10PercentFile() throws IOException{
        FileInputStream fin = new FileInputStream(new File("./src/test/resources/fastq/q30-10.fastq"));
        FastqReader fr = new FastqReader(fin);
        FastqData data = fr.read();
        double q30Percentage = data.getQ30() * 100.0 / data.getYield();
        assertEquals(100, data.getYield());
        assertEquals(10.0, q30Percentage, 0.00001);
    }
}

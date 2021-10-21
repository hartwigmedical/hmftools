package com.hartwig.hmftools.common.purple.cnchromosome;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class CnPerChromosomeArmFileTest {

    private static final String CNCHROMOSOME = Resources.getResource("purple/cnchromosome/sample.cnv.chromosomearm.somatic.tsv").getPath();


    @Test
    public void canReadFile() throws IOException {
        List<CnPerChromosomeArmData> cnPerChromosomeArmData =  CnPerChromosomeArmFile.read(CNCHROMOSOME);
        assertEquals(46, cnPerChromosomeArmData.size());
    }
}
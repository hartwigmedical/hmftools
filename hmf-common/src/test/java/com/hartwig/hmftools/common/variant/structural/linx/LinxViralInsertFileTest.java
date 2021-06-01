package com.hartwig.hmftools.common.variant.structural.linx;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class LinxViralInsertFileTest {

    private static final String LINX_VIRAL_INSERTIONS_FILE =
            Resources.getResource("variant/structural/linx/sample.linx.viral_inserts.tsv").getPath();

    @Test
    public void canReadViralInsertionFile() throws IOException {
        List<LinxViralInsertion> viralInsertionsFile = LinxViralInsertion.read(LINX_VIRAL_INSERTIONS_FILE);
        assertEquals(3, viralInsertionsFile.size());

        LinxViralInsertion viralInsertion1 = viralInsertionsFile.get(0);
        assertEquals("sample", viralInsertion1.SampleId);
        assertEquals(6, viralInsertion1.SvId);
        assertEquals("NC_001526", viralInsertion1.VirusId);
        assertEquals("Human papillomavirus type 16", viralInsertion1.VirusName);

        LinxViralInsertion viralInsertion2 = viralInsertionsFile.get(1);
        assertEquals("sample", viralInsertion2.SampleId);
        assertEquals(7, viralInsertion2.SvId);
        assertEquals("NC_001526", viralInsertion2.VirusId);
        assertEquals("Human papillomavirus type 16", viralInsertion2.VirusName);

        LinxViralInsertion viralInsertion3 = viralInsertionsFile.get(2);
        assertEquals("sample", viralInsertion3.SampleId);
        assertEquals(7, viralInsertion3.SvId);
        assertEquals("NC_001525", viralInsertion3.VirusId);
        assertEquals("Human papillomavirus type 15", viralInsertion3.VirusName);
    }
}
package com.hartwig.hmftools.patientreporter.viralInsertion;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.junit.Test;

public class InterpretViralInsertionTest {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String LINX_VIRAL_INSERTIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.viral_inserts.tsv";

    @Test
    public void canMergeViralInsertions() throws IOException {
        List<ViralInsertion> viralInsertions = InterpretViralInsertion.interpretVirals(LINX_VIRAL_INSERTIONS_TSV);
        assertEquals(2, viralInsertions.size());

        ViralInsertion viralInsertion1 = viralInsertions.get(0);
        assertEquals("Human papillomavirus type 15", viralInsertion1.virus());
        assertEquals(1, viralInsertion1.countVirus());

        ViralInsertion viralInsertion2 = viralInsertions.get(1);
        assertEquals("Human papillomavirus type 16", viralInsertion2.virus());
        assertEquals(2, viralInsertion2.countVirus());

    }

}
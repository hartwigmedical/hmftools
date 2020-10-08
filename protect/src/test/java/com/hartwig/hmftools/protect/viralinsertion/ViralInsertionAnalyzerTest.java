package com.hartwig.hmftools.protect.viralinsertion;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;

import org.junit.Test;

public class ViralInsertionAnalyzerTest {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String LINX_VIRAL_INSERTIONS_TSV = BASE_DIRECTORY + "/linx/sample.linx.viral_inserts.tsv";

    @Test
    public void canMergeViralInsertionsForReporting() throws IOException {
        List<LinxViralInsertion> linxViralInsertions = LinxViralInsertion.read(LINX_VIRAL_INSERTIONS_TSV);
        List<ViralInsertion> viralInsertions = ViralInsertionAnalyzer.analyzeViralInsertions(linxViralInsertions);
        assertEquals(2, viralInsertions.size());

        ViralInsertion viralInsertion1 = viralInsertions.get(0);
        assertEquals("Human papillomavirus type 15", viralInsertion1.virus());
        assertEquals(1, viralInsertion1.viralInsertionCount());

        ViralInsertion viralInsertion2 = viralInsertions.get(1);
        assertEquals("Human papillomavirus type 16", viralInsertion2.virus());
        assertEquals(2, viralInsertion2.viralInsertionCount());
    }
}
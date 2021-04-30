package com.hartwig.hmftools.protect.viralbreakend;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.protect.viralbreakend.ViralBreakendFactory;
import com.hartwig.hmftools.protect.viralbreakend.ViralInsertionAnalyzer;
import com.hartwig.hmftools.protect.viralbreakend.Viralbreakend;

import org.junit.Test;

public class ViralInsertionAnalyzerTest {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String VIRAL_BREAKEND_TSV = BASE_DIRECTORY + "/viralbreakend/sample_with.virusbreakend.vcf.summary.tsv";

    @Test
    public void canMergeViralInsertionsForReporting() throws IOException {
        List<Viralbreakend> linxViralInviralBreakend = ViralBreakendFactory.readViralBreakend(VIRAL_BREAKEND_TSV);
        List<Viralbreakend> viralBreakend = ViralInsertionAnalyzer.analyzeViralInsertions(linxViralInviralBreakend);
        assertEquals(1, viralBreakend.size());

        Viralbreakend viralBreakend1 = viralBreakend.get(0);
        assertEquals("Human papillomavirus type 16", viralBreakend1.nameAssigned());
        assertEquals("1", viralBreakend1.integrations());

    }
}
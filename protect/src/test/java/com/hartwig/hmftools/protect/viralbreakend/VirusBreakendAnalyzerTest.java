package com.hartwig.hmftools.protect.viralbreakend;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakend;
import com.hartwig.hmftools.common.virusbreakend.VirusBreakendFactory;

import org.junit.Test;

public class VirusBreakendAnalyzerTest {

    private static final String BASE_DIRECTORY = Resources.getResource("test_run").getPath();
    private static final String VIRUS_BREAKEND_TSV = BASE_DIRECTORY + "/virusbreakend/sample_with.virusbreakend.vcf.summary.tsv";

    @Test
    public void canMergeVirusBreakendsForReporting() throws IOException {
        List<VirusBreakend> virusBreakend = VirusBreakendFactory.readVirusBreakend(VIRUS_BREAKEND_TSV);
        List<VirusBreakend> virusBreakendAnalyzed = VirusBreakendAnalyzer.analyzeVirusBreakends(virusBreakend);
        assertEquals(1, virusBreakendAnalyzed.size());

        VirusBreakend viralBreakend1 = virusBreakendAnalyzed.get(0);
        assertEquals("Human papillomavirus type 16", viralBreakend1.nameAssigned());
        assertEquals(1, viralBreakend1.integrations());

    }
}
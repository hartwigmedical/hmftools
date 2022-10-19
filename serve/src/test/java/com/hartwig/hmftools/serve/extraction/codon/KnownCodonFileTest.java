package com.hartwig.hmftools.serve.extraction.codon;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class KnownCodonFileTest {

    @Test
    public void canReadFromFileAndConvert() throws IOException {
        String knownCodonTsv = KnownCodonFile.knownCodonTsvPath(ActionabilityTestUtil.TEST_SERVE_OUTPUT_DIR, RefGenomeVersion.V37);
        List<KnownCodon> knownCodons = KnownCodonFile.read(knownCodonTsv);

        assertEquals(2, knownCodons.size());

        List<String> lines = KnownCodonFile.toLines(knownCodons);
        List<KnownCodon> regeneratedCodons = KnownCodonFile.fromLines(lines);
        List<String> regeneratedLines = KnownCodonFile.toLines(regeneratedCodons);
        assertEquals(lines.size(), regeneratedLines.size());

        for (int i = 0; i < lines.size(); i++) {
            assertEquals(lines.get(i), regeneratedLines.get(i));
        }
    }
}
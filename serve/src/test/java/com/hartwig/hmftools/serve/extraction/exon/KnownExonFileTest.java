package com.hartwig.hmftools.serve.extraction.exon;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.serve.actionability.ActionabilityTestUtil;

import org.junit.Test;

public class KnownExonFileTest {

    @Test
    public void canReadFromFileAndConvert() throws IOException {
        String knownExonTsv = KnownExonFile.knownExonTsvPath(ActionabilityTestUtil.TEST_SERVE_OUTPUT_DIR, RefGenomeVersion.V37);
        List<KnownExon> knownExons = KnownExonFile.read(knownExonTsv);

        assertEquals(2, knownExons.size());

        List<String> lines = KnownExonFile.toLines(knownExons);
        List<KnownExon> regeneratedExons = KnownExonFile.fromLines(lines);
        List<String> regeneratedLines = KnownExonFile.toLines(regeneratedExons);
        assertEquals(lines.size(), regeneratedLines.size());

        for (int i = 0; i < lines.size(); i++) {
            assertEquals(lines.get(i), regeneratedLines.get(i));
        }
    }
}
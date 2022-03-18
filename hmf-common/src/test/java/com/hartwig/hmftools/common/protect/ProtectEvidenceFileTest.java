package com.hartwig.hmftools.common.protect;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ProtectEvidenceFileTest {

    private static final String EVIDENCE_TSV = Resources.getResource("protect/example.tsv").getPath();

    @Test
    public void canReadProtectEvidenceFile() throws IOException {
        List<ProtectEvidence> evidences = ProtectEvidenceFile.read(EVIDENCE_TSV);

        assertEquals(5, evidences.size());

        // Check one case with multiple sources and urls separately.
        ProtectEvidence evidence = findByTreatment(evidences, "Dabrafenib + Trametinib");
        assertEquals(2, evidence.sources().size());
        assertEquals(3, evidence.sourceUrls().size());
    }

    @NotNull
    private static ProtectEvidence findByTreatment(@NotNull List<ProtectEvidence> evidences, @NotNull String treatment) {
        for (ProtectEvidence evidence : evidences) {
            if (evidence.treatment().equals(treatment)) {
                return evidence;
            }
        }

        throw new IllegalStateException("Could not find evidence with treatment: " + treatment);
    }
}
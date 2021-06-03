package com.hartwig.hmftools.common.virus;

import static org.junit.Assert.assertEquals;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VirusBreakendQCStatusTest {

    @Test
    public void canExtractQCOfVirusBreakend() {
        assertEquals(VirusBreakendQCStatus.NO_ABNORMALITIES, VirusBreakendQCStatus.convert(Strings.EMPTY));
        assertEquals(VirusBreakendQCStatus.LOW_VIRAL_COVERAGE, VirusBreakendQCStatus.convert("LOW_VIRAL_COVERAGE"));
        assertEquals(VirusBreakendQCStatus.EXCESSIVE_VIRAL_COVERAGE,
                VirusBreakendQCStatus.convert("EXCESSIVE_VIRAL_COVERAGE"));
        assertEquals(VirusBreakendQCStatus.ASSEMBLY_DOWNSAMPLED,
                VirusBreakendQCStatus.convert("ASSEMBLY_DOWNSAMPLED"));
        assertEquals(VirusBreakendQCStatus.CHILD_TAXID_REFERENCE,
                VirusBreakendQCStatus.convert("CHILD_TAXID_REFERENCE"));
        assertEquals(VirusBreakendQCStatus.UNCLEAR_TAXID_ASSIGNMENT,
                VirusBreakendQCStatus.convert("UNCLEAR_TAXID_ASSIGNMENT"));
    }

    @Test(expected = IllegalStateException.class)
    public void crashOnUnknownQCStatusOfVirusBreakend() {
        VirusBreakendQCStatus.convert("ABC");
    }
}
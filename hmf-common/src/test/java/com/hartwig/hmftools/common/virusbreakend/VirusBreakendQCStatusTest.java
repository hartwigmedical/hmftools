package com.hartwig.hmftools.common.virusbreakend;

import static org.junit.Assert.assertEquals;

import org.apache.logging.log4j.util.Strings;
import org.junit.Test;

public class VirusBreakendQCStatusTest {

    @Test
    public void canExtractQCOfVirusBreakend() {
        assertEquals(VirusBreakendQCStatus.PASS, VirusBreakendQCStatus.extractVirusBreakendQCStatus(Strings.EMPTY));
        assertEquals(VirusBreakendQCStatus.LOW_VIRAL_COVERAGE, VirusBreakendQCStatus.extractVirusBreakendQCStatus("LOW_VIRAL_COVERAGE"));
        assertEquals(VirusBreakendQCStatus.EXCESSIVE_VIRAL_COVERAGE,
                VirusBreakendQCStatus.extractVirusBreakendQCStatus("EXCESSIVE_VIRAL_COVERAGE"));
        assertEquals(VirusBreakendQCStatus.ASSEMBLY_DOWNSAMPLED,
                VirusBreakendQCStatus.extractVirusBreakendQCStatus("ASSEMBLY_DOWNSAMPLED"));
        assertEquals(VirusBreakendQCStatus.CHILD_TAXID_REFERENCE,
                VirusBreakendQCStatus.extractVirusBreakendQCStatus("CHILD_TAXID_REFERENCE"));
        assertEquals(VirusBreakendQCStatus.UNCLEAR_TAXID_ASSIGNMENT,
                VirusBreakendQCStatus.extractVirusBreakendQCStatus("UNCLEAR_TAXID_ASSIGNMENT"));
    }

    @Test(expected = IllegalStateException.class)
    public void hasUnknownQCStatusOfVirusBreakend() {
        VirusBreakendQCStatus.extractVirusBreakendQCStatus("ABC");
    }
}
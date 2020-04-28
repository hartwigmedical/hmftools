package com.hartwig.hmftools.common.lims;

import static org.junit.Assert.*;

import org.junit.Test;

public class LimsViralInsertionChoiceTest {

    @Test
    public void canExtractLimsViralInsertionsChoice() {
        assertEquals(LimsViralInsertionChoice.REPORT_VIRAL_INSERION,
                LimsViralInsertionChoice.fromLimsViralInsertionsReportingChoiceString(true, "WIDE00000001T"));

        assertEquals(LimsViralInsertionChoice.REPORT_VIRAL_INSERION,
                LimsViralInsertionChoice.fromLimsViralInsertionsReportingChoiceString(false, "WIDE00000001T"));

        assertEquals(LimsViralInsertionChoice.NO_REPORT_VIRAL_INSERTIONS,
                LimsViralInsertionChoice.fromLimsViralInsertionsReportingChoiceString(false, "CPCT00000001T"));

        assertEquals(LimsViralInsertionChoice.NO_REPORT_VIRAL_INSERTIONS,
                LimsViralInsertionChoice.fromLimsViralInsertionsReportingChoiceString(true, "CPCT00000001T"));
    }

}
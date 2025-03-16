package com.hartwig.hmftools.pavereverse;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public final class DnaVariantsTest extends ReversePaveTestBase
{
    @Test
    public void calculateVariantsForSuppliedTranscriptId()
    {
        String transcriptId = "ENST00000322764";
        BaseSequenceChange bsc = reversePave.calculateDnaVariant("ZYX", transcriptId, "c.5C>A");
//        check(bsc, "A", "G", "7", 143381572 + 5); // coding start + 5
    }
}
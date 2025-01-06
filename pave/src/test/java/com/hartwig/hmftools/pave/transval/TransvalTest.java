package com.hartwig.hmftools.pave.transval;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import org.junit.Test;

public class TransvalTest extends TransvalTestBase
{
    @Test
    public void snv()
    {
        TransvalSNV record = transval.calculateSNV("MTOR:p.L2230V");

        assertEquals("ENST00000361445", record.TranscriptId);
        assertEquals("1", record.Chromosome);
//        assertEquals(11182158, record.Position);
        assertFalse(record.SpansMultipleExons);

//        TransvarSnvMnv snv = (TransvarSnvMnv) record.annotation();
//        assertEquals("A", snv.gdnaRef());
//        assertEquals("C", snv.gdnaAlt());
//        assertEquals("TTA", snv.referenceCodon());
//        assertEquals("GTA", snv.candidateCodons().get(0));
//        assertEquals("GTC", snv.candidateCodons().get(1));
//        assertEquals("GTG", snv.candidateCodons().get(2));
//        assertEquals("GTT", snv.candidateCodons().get(3));
    }
}

package com.hartwig.hmftools.telo;

import static com.hartwig.hmftools.telo.TeloConstants.CANONICAL_TELOMERE_SEQ;
import static com.hartwig.hmftools.telo.TeloUtils.hasTelomericContent;

import com.hartwig.hmftools.telo.TeloUtils;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import org.apache.commons.lang3.mutable.MutableBoolean;
import org.junit.Test;

public class TestTeloUtils
{
    @Test
    public void testTelomericContent()
    {
        String readBases = "AGCT" + CANONICAL_TELOMERE_SEQ + "AGCT" + CANONICAL_TELOMERE_SEQ + CANONICAL_TELOMERE_SEQ + "GG";
        assertTrue(hasTelomericContent(readBases));

        readBases = "AGCT" + CANONICAL_TELOMERE_SEQ + "AGCT" + CANONICAL_TELOMERE_SEQ + "GG";
        assertFalse(hasTelomericContent(readBases));
    }

    @Test
    public void testFullyTelomeric()
    {
        String readBasesG = "AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG";
        String readBasesC = "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC";

        assertTrue(TeloUtils.isFullyGTelomeric(readBasesG));
        assertFalse(TeloUtils.isFullyCTelomeric(readBasesG));

        assertFalse(TeloUtils.isFullyGTelomeric(readBasesC));
        assertTrue(TeloUtils.isFullyCTelomeric(readBasesC));

        //readBasesC = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCAAACCCAAACCCTAACCCAAACCCACACCCCCACACCAACCCCCACCCCCACCACAACACCCCCCCCCCCCCCCCCCCCACC";

    }



}

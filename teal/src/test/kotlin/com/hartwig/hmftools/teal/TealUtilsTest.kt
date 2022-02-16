package com.hartwig.hmftools.teal

import kotlin.test.*

class TealUtilsTest
{
    @Test
    fun testTelomericContent()
    {
        var readBases = "AGCT" + TealConstants.CANONICAL_TELOMERE_SEQ + "AGCT" + TealConstants.CANONICAL_TELOMERE_SEQ + TealConstants.CANONICAL_TELOMERE_SEQ + "GG"
        assertTrue(TealUtils.hasTelomericContent(readBases))
        readBases = "AGCT" + TealConstants.CANONICAL_TELOMERE_SEQ + "AGCT" + TealConstants.CANONICAL_TELOMERE_SEQ + "GG"
        assertFalse(TealUtils.hasTelomericContent(readBases))
    }

    @Test
    fun testLikelyTelomeric()
    {
        val readBasesG = "AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG"
        val readBasesC = "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC"
        assertTrue(TealUtils.isLikelyGTelomeric(readBasesG))
        assertFalse(TealUtils.isLikelyCTelomeric(readBasesG))
        assertFalse(TealUtils.isLikelyGTelomeric(readBasesC))
        assertTrue(TealUtils.isLikelyCTelomeric(readBasesC))

        //readBasesC = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCAAACCCAAACCCTAACCCAAACCCACACCCCCACACCAACCCCCACCCCCACCACAACACCCCCCCCCCCCCCCCCCCCACC";
    }
}
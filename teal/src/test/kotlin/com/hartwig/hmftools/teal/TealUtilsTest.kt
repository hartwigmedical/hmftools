package com.hartwig.hmftools.teal

import com.hartwig.hmftools.common.sequencing.SequencingType
import com.hartwig.hmftools.teal.tellength.FragmentClassifier
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
        val fragmentClassifier = FragmentClassifier(sequencingType = SequencingType.ILLUMINA)

        val readBasesG = "AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG"
        val readBasesC = "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC"
        assertTrue(fragmentClassifier.isLikelyGTelomeric(readBasesG))
        assertFalse(fragmentClassifier.isLikelyCTelomeric(readBasesG))
        assertFalse(fragmentClassifier.isLikelyGTelomeric(readBasesC))
        assertTrue(fragmentClassifier.isLikelyCTelomeric(readBasesC))

        //readBasesC = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCAAACCCAAACCCTAACCCAAACCCACACCCCCACACCAACCCCCACCCCCACCACAACACCCCCCCCCCCCCCCCCCCCACC";
    }

    @Test
    fun testLikelyTelomericSbx()
    {
        val fragmentClassifier = FragmentClassifier(sequencingType = SequencingType.SBX)

        val readBasesG = "AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG"
        val readBasesC = "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC"
        assertTrue(fragmentClassifier.isLikelyGTelomeric(readBasesG))
        assertFalse(fragmentClassifier.isLikelyCTelomeric(readBasesG))
        assertFalse(fragmentClassifier.isLikelyGTelomeric(readBasesC))
        assertTrue(fragmentClassifier.isLikelyCTelomeric(readBasesC))

        //readBasesC = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCAAACCCAAACCCTAACCCAAACCCACACCCCCACACCAACCCCCACCCCCACCACAACACCCCCCCCCCCCCCCCCCCCACC";
    }
}
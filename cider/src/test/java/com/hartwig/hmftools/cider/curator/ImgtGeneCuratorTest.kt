package com.hartwig.hmftools.cider.curator

import com.hartwig.hmftools.cider.blastn.BlastnMatch
import com.hartwig.hmftools.common.genome.region.Strand
import kotlin.test.*

class ImgtGeneCuratorTest
{
    @Test
    fun testCalcPositionMismatch()
    {
    }

    @Test
    fun testMatchToQueryGenomicLocationPosStrand()
    {
        var blastnMatch = BlastnMatch(
            querySeqLen = 30,
            subjectTitle = "Homo sapiens chromosome 12, GRCh38.p13 Primary Assembly",
            percentageIdent = 100.0,
            queryCoverage = 100.0,
            alignmentLength = 30,
            numMismatch = 0,
            numGapOpenings = 0,
            queryAlignStart = 1,
            queryAlignEnd = 30,
            subjectAlignStart = 1001,
            subjectAlignEnd = 1030,
            subjectFrame = Strand.FORWARD,
            expectedValue = 0.1,
            bitScore = 150.0,
            alignedPartOfQuerySeq = "",
            alignedPartOfSubjectSeq = ""
        )

        var genomicLocation = ImgtGeneCurator.matchToQueryGenomicLocation(blastnMatch)

        assertNotNull(genomicLocation)
        assertTrue(genomicLocation.inPrimaryAssembly)
        assertEquals("chr12", genomicLocation.chromosome)
        assertEquals(1001, genomicLocation.posStart)
        assertEquals(1030, genomicLocation.posEnd)
        assertEquals(Strand.FORWARD, genomicLocation.strand)

        // test align start / end changes
        blastnMatch = blastnMatch.copy(
            alignmentLength = 27,
            queryAlignStart = 2,
            queryAlignEnd = 28)

        genomicLocation = ImgtGeneCurator.matchToQueryGenomicLocation(blastnMatch)

        assertNotNull(genomicLocation)
        assertTrue(genomicLocation.inPrimaryAssembly)
        assertEquals("chr12", genomicLocation.chromosome)
        assertEquals(1000, genomicLocation.posStart)
        assertEquals(1032, genomicLocation.posEnd)
        assertEquals(Strand.FORWARD, genomicLocation.strand)
    }

    @Test
    fun testMatchToQueryGenomicLocationNegStrand()
    {
        var blastnMatch = BlastnMatch(
            querySeqLen = 30,
            subjectTitle = "Homo sapiens chromosome 12, GRCh38.p13 Primary Assembly",
            percentageIdent = 100.0,
            queryCoverage = 100.0,
            alignmentLength = 30,
            numMismatch = 0,
            numGapOpenings = 0,
            queryAlignStart = 1,
            queryAlignEnd = 30,
            subjectAlignStart = 1030,
            subjectAlignEnd = 1001,
            subjectFrame = Strand.REVERSE,
            expectedValue = 0.1,
            bitScore = 150.0,
            alignedPartOfQuerySeq = "",
            alignedPartOfSubjectSeq = ""
        )

        var genomicLocation = ImgtGeneCurator.matchToQueryGenomicLocation(blastnMatch)

        assertNotNull(genomicLocation)
        assertTrue(genomicLocation.inPrimaryAssembly)
        assertEquals("chr12", genomicLocation.chromosome)
        assertEquals(1001, genomicLocation.posStart)
        assertEquals(1030, genomicLocation.posEnd)
        assertEquals(Strand.REVERSE, genomicLocation.strand)

        // test align start / end changes
        blastnMatch = blastnMatch.copy(
            alignmentLength = 27,
            queryAlignStart = 2,
            queryAlignEnd = 28)

        genomicLocation = ImgtGeneCurator.matchToQueryGenomicLocation(blastnMatch)

        assertNotNull(genomicLocation)
        assertTrue(genomicLocation.inPrimaryAssembly)
        assertEquals("chr12", genomicLocation.chromosome)
        assertEquals(999, genomicLocation.posStart)
        assertEquals(1031, genomicLocation.posEnd)
        assertEquals(Strand.REVERSE, genomicLocation.strand)
    }
}
package com.hartwig.hmftools.cider.curator

import com.hartwig.hmftools.common.blastn.BlastnMatch
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
            30,
            "Homo sapiens chromosome 12, GRCh38.p13 Primary Assembly",
            100.0,
            100.0,
            30,
            0,
            0,
            1,
            30,
            1001,
            1030,
            Strand.FORWARD,
            0.1,
            150.0,
            "",
            ""
        )

        var genomicLocation = ImgtGeneCurator.matchToQueryGenomicLocation(blastnMatch)

        assertNotNull(genomicLocation)
        assertTrue(genomicLocation.inPrimaryAssembly)
        assertEquals("chr12", genomicLocation.chromosome)
        assertEquals(1001, genomicLocation.posStart)
        assertEquals(1030, genomicLocation.posEnd)
        assertEquals(Strand.FORWARD, genomicLocation.strand)

        // test align start / end changes
        blastnMatch = BlastnMatch(
            30,
            "Homo sapiens chromosome 12, GRCh38.p13 Primary Assembly",
            100.0,
            100.0,
            27,
            0,
            0,
            2,
            28,
            1001,
            1030,
            Strand.FORWARD,
            0.1,
            150.0,
            "",
            ""
        )

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
            30,
            "Homo sapiens chromosome 12, GRCh38.p13 Primary Assembly",
            100.0,
            100.0,
            30,
            0,
            0,
            1,
            30,
            1030,
            1001,
            Strand.REVERSE,
            0.1,
            150.0,
            "",
            ""
        )

        var genomicLocation = ImgtGeneCurator.matchToQueryGenomicLocation(blastnMatch)

        assertNotNull(genomicLocation)
        assertTrue(genomicLocation.inPrimaryAssembly)
        assertEquals("chr12", genomicLocation.chromosome)
        assertEquals(1001, genomicLocation.posStart)
        assertEquals(1030, genomicLocation.posEnd)
        assertEquals(Strand.REVERSE, genomicLocation.strand)

        // test align start / end changes
        blastnMatch = BlastnMatch(
            30,
            "Homo sapiens chromosome 12, GRCh38.p13 Primary Assembly",
            100.0,
            100.0,
            27,
            0,
            0,
            2,
            28,
            1030,
            1001,
            Strand.REVERSE,
            0.1,
            150.0,
            "",
            ""
        )

        genomicLocation = ImgtGeneCurator.matchToQueryGenomicLocation(blastnMatch)

        assertNotNull(genomicLocation)
        assertTrue(genomicLocation.inPrimaryAssembly)
        assertEquals("chr12", genomicLocation.chromosome)
        assertEquals(999, genomicLocation.posStart)
        assertEquals(1031, genomicLocation.posEnd)
        assertEquals(Strand.REVERSE, genomicLocation.strand)
    }
}
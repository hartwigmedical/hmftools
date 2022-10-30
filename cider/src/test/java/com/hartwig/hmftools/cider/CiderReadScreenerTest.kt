package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.genome.region.GenomeRegions
import com.hartwig.hmftools.common.genome.region.Strand
import htsjdk.samtools.SAMRecord
import junit.framework.TestCase
import org.junit.Test

class CiderReadScreenerTest
{
    companion object
    {
        const val MAX_FRAGMENT_LENGTH = 1000
    }

    @Test
    fun testExtrapolateReadOffsetAtRefPosition()
    {
        // 50 bases alignment
        val alignmentStart = 1001
        val alignmentEnd = 1050

        // create a SAM record with soft clip on both sides
        val record = SAMRecord(null)
        record.alignmentStart = alignmentStart
        // soft clips on both sides
        record.cigarString = "30S50M20S"
        // 100 bases here
        record.readString = "GACAACGCCAAGAACTCACTGTCTCTGCAAATGAATGACCTGCGAGTCGAAGACACGGCTGTGTATTACTGTGCGAGACCGAAATTTTATAGTAATGGCT"
        record.readNegativeStrandFlag = false
        //record.setBaseQualityString(qualities);
        record.mappingQuality = 20
        record.duplicateReadFlag = false
        record.readUnmappedFlag = false
        record.properPairFlag = true
        record.readPairedFlag = true

        // make sure we set it up correctly
        TestCase.assertEquals(31, record.getReadPositionAtReferencePosition(alignmentStart))
        TestCase.assertEquals(80, record.getReadPositionAtReferencePosition(alignmentEnd))
        var readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart)
        // start offset is 30, due to 30 left soft clip
        TestCase.assertEquals(30, readOffset)
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd)
        TestCase.assertEquals(79, readOffset)
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 30)
        // start offset is 0, due to 30 left soft clip
        TestCase.assertEquals(0, readOffset)
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd + 20)
        // start offset is 99, due to 20 right soft clip
        TestCase.assertEquals(99, readOffset)
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd + 30)
        // out of range
        TestCase.assertEquals(-1, readOffset)

        // test out of range
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 31)
        TestCase.assertEquals(-1, readOffset)

        // this is over the end, so should get -1
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd + 30)
        TestCase.assertEquals(-1, readOffset)
    }

    @Test
    fun testExtrapolateReadOffsetAtRefPositionSplice()
    {
        val alignmentStart = 100000
        val alignmentEnd = 108099

        // based on a RNA read we have come across
        val record = SAMRecord(null)
        record.referenceName = "X"
        record.alignmentStart = alignmentStart
        record.readString = "GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCA" +
                "GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCA"
        record.baseQualityString = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" +
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
        record.readNegativeStrandFlag = false
        record.mappingQuality = 20
        record.duplicateReadFlag = false
        record.readUnmappedFlag = false
        record.properPairFlag = true
        record.firstOfPairFlag = true
        record.readPairedFlag = true

        // following sets up multiple blocks
        // 1. 100000-10049 (50)
        // 2. 107050-107069 (20)
        // 3. 108070-108099 (30)
        record.cigarString = "20S50M7000N20M1000N30M80S"

        // there are multiple align blocks within
        TestCase.assertEquals(21, record.getReadPositionAtReferencePosition(alignmentStart))
        TestCase.assertEquals(120, record.getReadPositionAtReferencePosition(alignmentEnd))

        // try a point 30 bases before first alignment block
        var readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 30)
        TestCase.assertEquals(-1, readOffset) // 20S only so 30 bases is too many

        // try a point 20 bases before first alignment block, i.e. within the soft clipped region
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 20)
        TestCase.assertEquals(0, readOffset) // 20S - 20

        // try a point within the first match block of 50 bases
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, 100020)
        TestCase.assertEquals(40, readOffset) // 40 since it is 20 after the 20 soft clip

        // try a point 100 bases after first match block
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, 100150)
        TestCase.assertEquals(170, readOffset) // 20 soft clip + 150

        // 200 bases after first block, should get nothing since it will be outside range
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, 100250)
        TestCase.assertEquals(-1, readOffset)

        // match 50 bases before second block
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, 107000)
        TestCase.assertEquals(20, readOffset) // 20S + 50M - 50 bases

        // this will match 10 bases into second block
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, 107060)
        TestCase.assertEquals(80, readOffset) // 20S + 50M + 10 bases

        // this will match 30 after second block
        readOffset = CiderReadScreener.extrapolateReadOffsetAtRefPosition(record, 107100)
        TestCase.assertEquals(120, readOffset) // 20S + 50M + 20M + 30
    }

    @Test
    fun testIsRelevantToAnchorLocation1()
    {
        val readLength = 150
        val mappedLength = 100

        // positive strand
        val anchorRefStart = 10000
        val anchorRefEnd = 10030

        // this tests the function to work out if a read is potentially near and on the
        // correct side of the anchor genome location
        val anchorLocations = arrayOf(
            VJAnchorGenomeLocation(VJGeneType.IGHV, GenomeRegionStrand("1", anchorRefStart, anchorRefEnd, Strand.FORWARD)),
            VJAnchorGenomeLocation(VJGeneType.IGHJ, GenomeRegionStrand("1", anchorRefStart, anchorRefEnd, Strand.REVERSE))
        )
        for (anchorLocation in anchorLocations)
        {
            // we should only allow reads that are mapped at lower genome location, i.e.
            // read ------ anchor ----CDR3
            // or reads that overlap the anchor by at least 15 bases

            // first try reads that are lower
            var mappedEnd = anchorRefEnd - readLength + 15
            var mappedStart = mappedEnd - mappedLength
            var mapped = GenomeRegions.create(anchorLocation.chromosome, mappedStart, mappedEnd)
            TestCase.assertTrue(CiderReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation))

            // reads with coords above and not overlapping anchor are not relevant
            mappedStart = anchorRefEnd + anchorLocation.baseLength()
            mappedEnd = mappedStart + mappedLength
            mapped = GenomeRegions.create(anchorLocation.chromosome, mappedStart, mappedEnd)
            TestCase.assertFalse(CiderReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation))

            // reads that overlap with anchor by 15 bases or more are relevant
            mappedEnd = anchorRefEnd + anchorLocation.baseLength() / 2
            mappedStart = mappedEnd - mappedLength
            mapped = GenomeRegions.create(anchorLocation.chromosome, mappedStart, mappedEnd)
            TestCase.assertTrue(CiderReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation))
        }
    }

    @Test
    fun testIsRelevantToAnchorLocation2()
    {
        val readLength = 150
        val mappedLength = 100

        // positive strand
        val anchorRefStart = 10000
        val anchorRefEnd = 10030

        // this tests the function to work out if a read is potentially near and on the
        // correct side of the anchor genome location
        val anchorLocations = arrayOf(
            VJAnchorGenomeLocation(VJGeneType.TRAV, GenomeRegionStrand("1", anchorRefStart, anchorRefEnd, Strand.REVERSE)),
            VJAnchorGenomeLocation(VJGeneType.TRAJ, GenomeRegionStrand("1", anchorRefStart, anchorRefEnd, Strand.FORWARD))
        )
        for (anchorLocation in anchorLocations)
        {
            // we should only allow reads that are mapped at higher genome location, i.e.
            // CDR3 ------ anchor ----read
            // or reads that overlap the anchor by at least 15 bases

            // first try reads that are mapped at higher coord
            var mappedStart = anchorRefStart + readLength - 15
            var mappedEnd = mappedStart + mappedLength
            var mapped = GenomeRegions.create(anchorLocation.chromosome, mappedStart, mappedEnd)
            TestCase.assertTrue(CiderReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation))

            // reads with coords below and not overlapping anchor are not relevant
            mappedEnd = anchorRefStart - anchorLocation.baseLength()
            mappedStart = mappedEnd - mappedLength
            mapped = GenomeRegions.create(anchorLocation.chromosome, mappedStart, mappedEnd)
            TestCase.assertFalse(CiderReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation))

            // reads that overlap with anchor by 15 bases or more are relevant
            mappedStart = anchorRefStart - anchorLocation.baseLength() / 2
            mappedEnd = mappedStart + mappedLength
            mapped = GenomeRegions.create(anchorLocation.chromosome, mappedStart, mappedEnd)
            TestCase.assertTrue(CiderReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation))
        }
    }

    @Test
    fun testFindAnchorPositionRNA()
    {
        // read: A00624:61:HVW7TDSXX:4:2344:18674:2018 1/2 151b aligned to 14:106322274-106330832., aligned: 121
        // A00624:61:HVW7TDSXX:4:2344:18674:2018	99	14	106322274	255	49M7085N20M1389N16M66S	=	106518456	196235
        // GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCAAGCCATTACTATAAAATTTCGGTCTCGCACAGTAATACACAGCCGTGTCTT
        // FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        // NH:i:1	HI:i:1	AS:i:120	nM:i:3	NM:i:2	MD:Z:50T2G31	jM:B:c,2,0	jI:B:i,106322323,106329407,106329428,106330816	MC:Z:16S53M82S

        // 50 bases alignment
        val alignmentStart = 106322274
        val alignmentEnd = 106330832

        // based on a RNA read we have come across
        val record = SAMRecord(null)
        record.readName = "read1"
        record.referenceName = "14"
        record.alignmentStart = alignmentStart
        record.cigarString = "49M7085N20M1389N16M66S"
        record.readString =
            "GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCAAGCCATTACTATAAAATTTCGGTCTCGCACAGTAATACACAGCCGTGTCTT"
        record.baseQualityString =
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
        record.readNegativeStrandFlag = false
        record.mappingQuality = 20
        record.duplicateReadFlag = false
        record.readUnmappedFlag = false
        record.properPairFlag = true
        record.firstOfPairFlag = true
        record.readPairedFlag = true
        val ighJ1 = VJAnchorTemplate(
            VJGeneType.IGHJ,
            "IGHJ1",
            "01",
            GenomeRegionStrand("14", 106330701, 106330840, Strand.REVERSE),
            "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG",
            "TGGGGCCAGGGCACCCTGGTCACCGTCTCC",
            GenomeRegionStrand("14", 106330801, 106330830, Strand.REVERSE)
        )
        val vjGeneStore = TestCiderGeneDatastore(listOf(ighJ1))
        val mockAnchorBlosumSearcher = MockAnchorBlosumSearcher()
        val ciderReadScreener = CiderReadScreener(
            vjGeneStore, mockAnchorBlosumSearcher, 6, MAX_FRAGMENT_LENGTH
        )

        // make sure we set it up correctly
        TestCase.assertEquals(1, record.getReadPositionAtReferencePosition(alignmentStart))
        TestCase.assertEquals(85, record.getReadPositionAtReferencePosition(alignmentEnd))
        val mapped = if (record.readUnmappedFlag) null else GenomeRegions.create(
            record.referenceName,
            record.alignmentStart, record.alignmentEnd
        )

        // template loc: 14:106330801-106330830(-)
        val anchorLocation = VJAnchorGenomeLocation(
            VJGeneType.IGHJ,
            GenomeRegionStrand("14", 106330801, 106330830, Strand.REVERSE)
        )
        val readCandidate = ciderReadScreener.matchesAnchorLocation(record, mapped!!, anchorLocation, true)
        TestCase.assertNotNull(readCandidate)
        TestCase.assertEquals(68, readCandidate!!.anchorOffsetStart)
        TestCase.assertEquals(98, readCandidate.anchorOffsetEnd)
    }

    // we test unmapped read where the mate is mapped to upstream of V
    @Test
    fun isUnamppedReadRelevantToAnchorLocV()
    {
        val chr = "1"
        // positive strand
        val anchorRefStart = 10000
        val anchorRefEnd = 10030

        // create a SAM record with soft clip on both sides
        val read = SAMRecord(null)
        read.mateReferenceName = chr
        read.cigarString = "*"
        // 100 bases here
        read.readString = "GACAACGCCAAGAACTCACTGTCTCTGCAAATGAATGACCTGCGAGTCGAAGACACGGCTGTGTATTACTGTGCGAGACCGAAATTTTATAGTAATGGCT"
        read.readNegativeStrandFlag = false
        read.readUnmappedFlag = true
        read.mateUnmappedFlag = false
        read.properPairFlag = true
        read.readPairedFlag = true

        // first test forward strand

        // Here is what we want to test:
        // 1. the mapped mate must have coord below the V anchor coord, but not too far away (max 1000 bases)
        // 2. the anchor location is forward strand
        // 3. the mapped mate must be mapped to forward strand
        // >------------V--------D---------J---------> forward strand
        // ======>  <=====
        //  mate     this
        //
        var anchorLocation = VJAnchorGenomeLocation(VJGeneType.TRAV, GenomeRegionStrand(chr, anchorRefStart, anchorRefEnd, Strand.FORWARD))

        read.mateAlignmentStart = anchorRefStart - 500
        read.mateNegativeStrandFlag = false
        TestCase.assertTrue(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // too far away
        read.mateAlignmentStart = anchorRefStart - MAX_FRAGMENT_LENGTH - read.readLength
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // if the mate is mapped to the other side of V, we do not accept
        read.mateAlignmentStart = anchorRefEnd + 500
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // if mate is negative strand, it is pointing at wrong direction so also do not take it
        read.mateAlignmentStart = anchorRefStart - 500
        read.mateNegativeStrandFlag = true
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // now test reverse strand
        // Here is what we want to test:
        // 1. the mapped mate must have coord higher the V anchor coord, but not too far away (max 1000 bases)
        // 2. the anchor location is reverse strand
        // 3. the mapped mate must be mapped to reverse strand
        // <-------J--------D--------V-------------< reverse strand
        //                          ======>  <=====
        //                           this     mate
        //
        anchorLocation = VJAnchorGenomeLocation(VJGeneType.TRAV, GenomeRegionStrand(chr, anchorRefStart, anchorRefEnd, Strand.REVERSE))

        read.mateAlignmentStart = anchorRefEnd + 500
        read.mateNegativeStrandFlag = true
        TestCase.assertTrue(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // too far away
        read.mateAlignmentStart = anchorRefEnd + MAX_FRAGMENT_LENGTH
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // if the mate is mapped to the other side of V, we do not accept
        read.mateAlignmentStart = anchorRefStart - 500
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // if mate is negative strand, it is pointing at wrong direction so also do not take it
        read.mateAlignmentStart = anchorRefEnd + 500
        read.mateNegativeStrandFlag = false
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))
    }

    // we test unmapped read where the mate is mapped to downstream of J
    @Test
    fun isUnamppedReadRelevantToAnchorLocJ()
    {
        val chr = "1"
        // positive strand
        val anchorRefStart = 10000
        val anchorRefEnd = 10030

        // create a SAM record with soft clip on both sides
        val read = SAMRecord(null)
        read.mateReferenceName = chr
        read.cigarString = "*"
        // 100 bases here
        read.readString = "GACAACGCCAAGAACTCACTGTCTCTGCAAATGAATGACCTGCGAGTCGAAGACACGGCTGTGTATTACTGTGCGAGACCGAAATTTTATAGTAATGGCT"
        read.readNegativeStrandFlag = false
        read.readUnmappedFlag = true
        read.mateUnmappedFlag = false
        read.properPairFlag = true
        read.readPairedFlag = true

        // first test forward strand

        // Here is what we want to test:
        // 1. the mapped mate must have coord above the J anchor coord, but not too far away (max 1000 bases)
        // 2. the anchor location is forward strand
        // 3. the mapped mate must be mapped to negative strand
        // >----V-------D---------J---------> forward strand
        //                      ======>  <=====
        //                       this     mate
        //
        var anchorLocation = VJAnchorGenomeLocation(VJGeneType.TRAJ, GenomeRegionStrand(chr, anchorRefStart, anchorRefEnd, Strand.FORWARD))

        read.mateAlignmentStart = anchorRefEnd + 500
        read.mateNegativeStrandFlag = true
        TestCase.assertTrue(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // too far away
        read.mateAlignmentStart = anchorRefEnd + MAX_FRAGMENT_LENGTH
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // if the mate is mapped to the other side of J, we do not accept
        read.mateAlignmentStart = anchorRefStart - 500
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // if mate is positive strand, it is pointing at wrong direction so also do not take it
        read.mateAlignmentStart = anchorRefEnd + 500
        read.mateNegativeStrandFlag = false
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // now test reverse strand
        // Here is what we want to test:
        // 1. the mapped mate must have coord higher the V anchor coord, but not too far away (max 1000 bases)
        // 2. the anchor location is reverse strand
        // 3. the mapped mate must be mapped to positive strand
        // <-----------J-------D----------V-------------< reverse strand
        // ======>  <=====
        //  mate     this
        //
        anchorLocation = VJAnchorGenomeLocation(VJGeneType.TRAJ, GenomeRegionStrand(chr, anchorRefStart, anchorRefEnd, Strand.REVERSE))

        read.mateAlignmentStart = anchorRefStart - 500
        read.mateNegativeStrandFlag = false
        TestCase.assertTrue(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // too far away
        read.mateAlignmentStart = anchorRefStart + MAX_FRAGMENT_LENGTH - read.readLength
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // if the mate is mapped to the other side of J, we do not accept
        read.mateAlignmentStart = anchorRefEnd + 500
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))

        // if mate is negative strand, it is pointing at wrong direction so also do not take it
        read.mateAlignmentStart = anchorRefStart - 500
        read.mateNegativeStrandFlag = true
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToAnchorLoc(read, anchorLocation, MAX_FRAGMENT_LENGTH))
    }

    // we test unmapped read where the mate is mapped around constant region
    @Test
    fun isUnamppedReadRelevantToConstantRegion()
    {
        val chr = "1"
        // positive strand
        val constantRegionRefStart = 10000
        val constantRegionRefEnd = 10030

        // create a SAM record with soft clip on both sides
        val read = SAMRecord(null)
        read.mateReferenceName = chr
        read.cigarString = "*"
        // 100 bases here
        read.readString = "GACAACGCCAAGAACTCACTGTCTCTGCAAATGAATGACCTGCGAGTCGAAGACACGGCTGTGTATTACTGTGCGAGACCGAAATTTTATAGTAATGGCT"
        read.readNegativeStrandFlag = false
        read.readUnmappedFlag = true
        read.mateUnmappedFlag = false
        read.properPairFlag = true
        read.readPairedFlag = true

        // first test forward strand

        // Here is what we want to test:
        // 1. the mapped mate must have coord near constant region
        // 2. the constant region is on forward strand
        // 3. the mapped mate must be mapped to negative strand
        // >----V---D---J------C----> forward strand
        //         ======>  <=====
        //          this     mate
        //
        var constantRegion = GenomeRegionStrand(chr, constantRegionRefStart, constantRegionRefEnd, Strand.FORWARD)

        read.mateAlignmentStart = constantRegionRefEnd + 500
        read.mateNegativeStrandFlag = true
        TestCase.assertTrue(CiderReadScreener.isUnamppedReadRelevantToConstantRegion(read, constantRegion, MAX_FRAGMENT_LENGTH))

        // the other side is also fine, we do not impose any rules for now
        read.mateAlignmentStart = constantRegionRefStart - 500
        TestCase.assertTrue(CiderReadScreener.isUnamppedReadRelevantToConstantRegion(read, constantRegion, MAX_FRAGMENT_LENGTH))

        // too far away
        read.mateAlignmentStart = constantRegionRefEnd + MAX_FRAGMENT_LENGTH
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToConstantRegion(read, constantRegion, MAX_FRAGMENT_LENGTH))

        // if mate is positive strand, it is pointing at wrong direction so also do not take it
        read.mateAlignmentStart = constantRegionRefEnd + 500
        read.mateNegativeStrandFlag = false
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToConstantRegion(read, constantRegion, MAX_FRAGMENT_LENGTH))

        // now test reverse strand

        // Here is what we want to test:
        // 1. the mapped mate must have coord near constant region
        // 2. the constant region is on reverse strand
        // 3. the mapped mate must be mapped to positive strand
        // <---C-------J---D---V------< reverse strand
        //   ======>   <=====
        //     mate     this
        //
        constantRegion = GenomeRegionStrand(chr, constantRegionRefStart, constantRegionRefEnd, Strand.REVERSE)

        read.mateAlignmentStart = constantRegionRefStart - 500
        read.mateNegativeStrandFlag = false
        TestCase.assertTrue(CiderReadScreener.isUnamppedReadRelevantToConstantRegion(read, constantRegion, MAX_FRAGMENT_LENGTH))

        // the other side is also fine, we do not impose any rule for now
        read.mateAlignmentStart = constantRegionRefEnd + 500
        TestCase.assertTrue(CiderReadScreener.isUnamppedReadRelevantToConstantRegion(read, constantRegion, MAX_FRAGMENT_LENGTH))

        // too far away
        read.mateAlignmentStart = constantRegionRefStart - MAX_FRAGMENT_LENGTH - read.readLength
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToConstantRegion(read, constantRegion, MAX_FRAGMENT_LENGTH))

        // if mate is negative strand, it is pointing at wrong direction so also do not take it
        read.mateAlignmentStart = constantRegionRefStart - 500
        read.mateNegativeStrandFlag = true
        TestCase.assertFalse(CiderReadScreener.isUnamppedReadRelevantToConstantRegion(read, constantRegion, MAX_FRAGMENT_LENGTH))
    }
}
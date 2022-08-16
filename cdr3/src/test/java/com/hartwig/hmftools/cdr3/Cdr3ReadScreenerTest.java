package com.hartwig.hmftools.cdr3;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import java.util.HashSet;
import java.util.List;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class Cdr3ReadScreenerTest
{
    @Test
    public void testExtrapolateReadOffsetAtRefPosition()
    {
        // 50 bases alignment
        int alignmentStart = 1001;
        int alignmentEnd = 1050;

        // create a SAM record with soft clip on both sides
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(alignmentStart);
        // soft clips on both sides
        record.setCigarString("30S50M20S");
        // 100 bases here
        record.setReadString("GACAACGCCAAGAACTCACTGTCTCTGCAAATGAATGACCTGCGAGTCGAAGACACGGCTGTGTATTACTGTGCGAGACCGAAATTTTATAGTAATGGCT");
        record.setReadNegativeStrandFlag(false);
        //record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setReadPairedFlag(true);

        // make sure we set it up correctly
        assertEquals(31, record.getReadPositionAtReferencePosition(alignmentStart));
        assertEquals(80, record.getReadPositionAtReferencePosition(alignmentEnd));

        int readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart);
        // start offset is 30, due to 30 left soft clip
        assertEquals(30, readOffset);

        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd);
        assertEquals(79, readOffset);

        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 30);
        // start offset is 0, due to 30 left soft clip
        assertEquals(0, readOffset);

        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd + 20);
        // start offset is 99, due to 20 right soft clip
        assertEquals(99, readOffset);

        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd + 30);
        // out of range
        assertEquals(-1, readOffset);

        // test out of range
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 31);
        assertEquals(-1, readOffset);

        // this is over the end, so should get -1
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentEnd + 30);
        assertEquals(-1, readOffset);
    }

    @Test
    public void testExtrapolateReadOffsetAtRefPositionSplice()
    {
        int alignmentStart = 100000;
        int alignmentEnd   = 108099;

        // based on a RNA read we have come across
        final SAMRecord record = new SAMRecord(null);
        record.setReferenceName("X");
        record.setAlignmentStart(alignmentStart);
        record.setReadString("GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCA" +
                "GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCA");
        record.setBaseQualityString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" +
                "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
        record.setReadNegativeStrandFlag(false);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setFirstOfPairFlag(true);
        record.setReadPairedFlag(true);

        // following sets up multiple blocks
        // 1. 100000-10049 (50)
        // 2. 107050-107069 (20)
        // 3. 108070-108099 (30)
        record.setCigarString("20S50M7000N20M1000N30M80S");

        // there are multiple align blocks within
        assertEquals(21, record.getReadPositionAtReferencePosition(alignmentStart));
        assertEquals(120, record.getReadPositionAtReferencePosition(alignmentEnd));

        // try a point 30 bases before first alignment block
        int readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 30);
        assertEquals(-1, readOffset); // 20S only so 30 bases is too many

        // try a point 20 bases before first alignment block, i.e. within the soft clipped region
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, alignmentStart - 20);
        assertEquals(0, readOffset); // 20S - 20

        // try a point within the first match block of 50 bases
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, 100020);
        assertEquals(40, readOffset); // 40 since it is 20 after the 20 soft clip

        // try a point 100 bases after first match block
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, 100150);
        assertEquals(170, readOffset); // 20 soft clip + 150

        // 200 bases after first block, should get nothing since it will be outside range
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, 100250);
        assertEquals(-1, readOffset);

        // match 50 bases before second block
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, 107000);
        assertEquals(20, readOffset); // 20S + 50M - 50 bases

        // this will match 10 bases into second block
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, 107060);
        assertEquals(80, readOffset); // 20S + 50M + 10 bases

        // this will match 30 after second block
        readOffset = Cdr3ReadScreener.extrapolateReadOffsetAtRefPosition(record, 107100);
        assertEquals(120, readOffset); // 20S + 50M + 20M + 30
    }

    @Test
    public void testIsRelevantToAnchorLocation1()
    {
        int readLength = 150;
        int mappedLength = 100;

        // positive strand
        int anchorRefStart = 10000;
        int anchorRefEnd = 10030;

        // this tests the function to work out if a read is potentially near and on the
        // correct side of the anchor genome location
        VJAnchorReferenceLocation[] anchorLocations = {
                new VJAnchorReferenceLocation(VJ.V, new GeneLocation("1", anchorRefStart, anchorRefEnd, Strand.FORWARD)),
                new VJAnchorReferenceLocation(VJ.J, new GeneLocation("1", anchorRefStart, anchorRefEnd, Strand.REVERSE)),
        };

        for (VJAnchorReferenceLocation anchorLocation : anchorLocations)
        {
            // we should only allow reads that are mapped at lower genome location, i.e.
            // read ------ anchor ----CDR3
            // or reads that overlap the anchor by at least 15 bases

            // first try reads that are lower
            int mappedEnd = anchorRefEnd - readLength + 15;
            int mappedStart = mappedEnd - mappedLength;
            GenomeRegion mapped = GenomeRegions.create(anchorLocation.getChromosome(), mappedStart, mappedEnd);
            assertTrue(Cdr3ReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation));

            // reads with coords above and not overlapping anchor are not relevant
            mappedStart = anchorRefEnd + anchorLocation.baseLength();
            mappedEnd = mappedStart + mappedLength;
            mapped = GenomeRegions.create(anchorLocation.getChromosome(), mappedStart, mappedEnd);
            assertFalse(Cdr3ReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation));

            // reads that overlap with anchor by 15 bases or more are relevant
            mappedEnd = anchorRefEnd + anchorLocation.baseLength() / 2;
            mappedStart = mappedEnd - mappedLength;
            mapped = GenomeRegions.create(anchorLocation.getChromosome(), mappedStart, mappedEnd);
            assertTrue(Cdr3ReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation));
        }
    }

    @Test
    public void testIsRelevantToAnchorLocation2()
    {
        int readLength = 150;
        int mappedLength = 100;

        // positive strand
        int anchorRefStart = 10000;
        int anchorRefEnd = 10030;

        // this tests the function to work out if a read is potentially near and on the
        // correct side of the anchor genome location
        VJAnchorReferenceLocation[] anchorLocations = {
                new VJAnchorReferenceLocation(VJ.V, new GeneLocation("1", anchorRefStart, anchorRefEnd, Strand.REVERSE)),
                new VJAnchorReferenceLocation(VJ.J, new GeneLocation("1", anchorRefStart, anchorRefEnd, Strand.FORWARD)),
        };

        for (VJAnchorReferenceLocation anchorLocation : anchorLocations)
        {
            // we should only allow reads that are mapped at higher genome location, i.e.
            // CDR3 ------ anchor ----read
            // or reads that overlap the anchor by at least 15 bases

            // first try reads that are mapped at higher coord
            int mappedStart = anchorRefStart + readLength - 15;
            int mappedEnd = mappedStart + mappedLength;
            GenomeRegion mapped = GenomeRegions.create(anchorLocation.getChromosome(), mappedStart, mappedEnd);
            assertTrue(Cdr3ReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation));

            // reads with coords below and not overlapping anchor are not relevant
            mappedEnd = anchorRefStart - anchorLocation.baseLength();
            mappedStart = mappedEnd - mappedLength;
            mapped = GenomeRegions.create(anchorLocation.getChromosome(), mappedStart, mappedEnd);
            assertFalse(Cdr3ReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation));

            // reads that overlap with anchor by 15 bases or more are relevant
            mappedStart = anchorRefStart - anchorLocation.baseLength() / 2;
            mappedEnd = mappedStart + mappedLength;
            mapped = GenomeRegions.create(anchorLocation.getChromosome(), mappedStart, mappedEnd);
            assertTrue(Cdr3ReadScreener.isRelevantToAnchorLocation(readLength, mapped, anchorLocation));
        }
    }

    @Test
    public void testFindAnchorPositionRNA()
    {
        // read: A00624:61:HVW7TDSXX:4:2344:18674:2018 1/2 151b aligned to 14:106322274-106330832., aligned: 121
        // A00624:61:HVW7TDSXX:4:2344:18674:2018	99	14	106322274	255	49M7085N20M1389N16M66S	=	106518456	196235
        // GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCAAGCCATTACTATAAAATTTCGGTCTCGCACAGTAATACACAGCCGTGTCTT
        // FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        // NH:i:1	HI:i:1	AS:i:120	nM:i:3	NM:i:2	MD:Z:50T2G31	jM:B:c,2,0	jI:B:i,106322323,106329407,106329428,106330816	MC:Z:16S53M82S

        // 50 bases alignment
        int alignmentStart = 106322274;
        int alignmentEnd = 106330832;

        // based on a RNA read we have come across
        final SAMRecord record = new SAMRecord(null);
        record.setReadName("read1");
        record.setReferenceName("14");
        record.setAlignmentStart(alignmentStart);
        record.setCigarString("49M7085N20M1389N16M66S");
        record.setReadString("GAATTCTCACAGGAGACGAGGGGGAAAAGGGTTGGGGCGGATGCACTCCCGGAAGAGACGGTGACCGTGGTCCCTTGGCCCCAGACGTCCATACCCGCCAAGCCATTACTATAAAATTTCGGTCTCGCACAGTAATACACAGCCGTGTCTT");
        record.setBaseQualityString("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF");
        record.setReadNegativeStrandFlag(false);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        record.setProperPairFlag(true);
        record.setFirstOfPairFlag(true);
        record.setReadPairedFlag(true);

        var ighJ1 = new VJGene(
                "IGHJ1*01",
                "IGHJ1",
                "01",
                new GeneLocation("14", 106330701, 106330840, Strand.REVERSE),
                "GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG",
                "TGGGGCCAGGGCACCCTGGTCACCGTCTCC",
                new GeneLocation("14", 106330801, 106330830, Strand.REVERSE));
        var vjGeneStore = new TestVJGeneStore(List.of(ighJ1));
        var mockAnchorBlosumSearcher = new MockAnchorBlosumSearcher();

        Cdr3ReadScreener cdr3ReadScreener = new Cdr3ReadScreener(vjGeneStore, mockAnchorBlosumSearcher,
                151, 6);

        // make sure we set it up correctly
        assertEquals(1, record.getReadPositionAtReferencePosition(alignmentStart));
        assertEquals(85, record.getReadPositionAtReferencePosition(alignmentEnd));

        GenomeRegion mapped = record.getReadUnmappedFlag() ? null : GenomeRegions.create(record.getReferenceName(),
                record.getAlignmentStart(), record.getAlignmentEnd());

        // template loc: 14:106330801-106330830(-)
        VJAnchorReferenceLocation anchorLocation = new VJAnchorReferenceLocation(
                VJ.J,
                new GeneLocation("14", 106330801, 106330830, Strand.REVERSE));

        VJReadCandidate readCandidate = cdr3ReadScreener.matchesAnchorLocation(record, mapped, new HashSet<>(), anchorLocation);

        assertNotNull(readCandidate);
        assertEquals(68, readCandidate.getAnchorOffsetStart());
        assertEquals(98, readCandidate.getAnchorOffsetEnd());
    }
}

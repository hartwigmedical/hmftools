package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.buildDefaultBaseQuals;
import static com.hartwig.hmftools.sage.common.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.sage.common.TestUtils.buildCigarString;
import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadContext;
import static com.hartwig.hmftools.sage.common.VariantUtils.createReadCounter;
import static com.hartwig.hmftools.sage.evidence.SplitReadSegment.formSegment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.sage.common.VariantReadContext;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SplitReadTest
{
    @Test
    public void testSplitReadSegments()
    {
        String readBases = REF_BASES_200.substring(0, 50);
        String readCigar = "20M100N30M";
        byte[] readQuals = buildDefaultBaseQuals(readBases.length());

        // 90-109 = 20M, 110 - 209 = 100D, 210 - 239 = 30M
        SAMRecord read = buildSamRecord(90, readCigar, readBases, readQuals);

        SplitReadSegment splitReadSegment = formSegment(read, 100, 20);

        assertNotNull(splitReadSegment);
        assertEquals(20, splitReadSegment.ReadVarIndex);
        assertEquals(0, splitReadSegment.SegmentIndexStart);
        assertEquals(19, splitReadSegment.SegmentIndexEnd);
        assertEquals(20, splitReadSegment.ReadBases.length);

        splitReadSegment = formSegment(read, 115, 20);

        assertNull(splitReadSegment);

        splitReadSegment = formSegment(read, 215, 26);

        assertNotNull(splitReadSegment);
        assertEquals(6, splitReadSegment.ReadVarIndex);
        assertEquals(20, splitReadSegment.SegmentIndexStart);
        assertEquals(49, splitReadSegment.SegmentIndexEnd);
        assertEquals(30, splitReadSegment.ReadBases.length);

        // with a second split and soft-clips at the start and end
        readBases = REF_BASES_200.substring(0, 70);
        readCigar = "5S20M100N20M100N20M5S";
        readQuals = buildDefaultBaseQuals(readBases.length());

        // 85-89 = 5S, 90-109 = 20M, 110-209 = 100D, 210-229 = 20M, 230-329 = 100N, 330-349 = 20M, 350-354 = 5S
        read = buildSamRecord(90, readCigar, readBases, readQuals);

        splitReadSegment = formSegment(read, 335, 51);

        assertNotNull(splitReadSegment);
        assertEquals(6, splitReadSegment.ReadVarIndex);
        assertEquals(45, splitReadSegment.SegmentIndexStart);
        assertEquals(69, splitReadSegment.SegmentIndexEnd);
        assertEquals(25, splitReadSegment.ReadBases.length);

        // with normal indels as well
        readBases = REF_BASES_200.substring(0, 60);
        readCigar = "10M100N10M10D10M10I10M100N10M";
        readQuals = buildDefaultBaseQuals(readBases.length());

        // 90-99 = 10M, 100-199 = 100D, 200-209 = 10M, 210-219 = 10D, 220-229 = 10M, 229 = 10I, 230-239 = 10M, 240-339 = 100N, 340-349 = 10M
        read = buildSamRecord(90, readCigar, readBases, readQuals);

        splitReadSegment = formSegment(read, 235, 46);

        assertNotNull(splitReadSegment);
        assertEquals(36, splitReadSegment.ReadVarIndex);
        assertEquals(10, splitReadSegment.SegmentIndexStart);
        assertEquals(49, splitReadSegment.SegmentIndexEnd);
        assertEquals(40, splitReadSegment.ReadBases.length);
    }

    @Test
    public void testSplitReads()
    {
        int position = 133;
        VariantReadContext readContext = createReadContext(position, "A", "T");

        ReadContextCounter readContextCounter = createReadCounter(0, readContext);

        // first a read after an N section
        String altReadBases = REF_BASES_200.substring(0, 10) + readContext.readBases() + REF_BASES_200.substring(0, 10);
        assertEquals(45, altReadBases.length());

        String readBases = REF_BASES_200.substring(0, 10) + altReadBases + REF_BASES_200.substring(0, 10);

        // 1-10 = 10M, 11-110 = 100N, 111-155 = 45M, 156-255 = 100N, 256-265 = 10M
        // with variant at 210-219 buffer, 220-229 left flank, 230-234 core, so position at 232, read index 32
        String readCigar = "10M100N45M100N10M";

        SAMRecord read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 1, readBases, readCigar);

        readContextCounter.processRead(read, 1, null);
        assertEquals(1, readContextCounter.readCounts().Full);

        // now a read just covering the core
        readCigar = "10M122N23M100N10M";

        readBases = REF_BASES_200.substring(0, 10) + altReadBases.substring(22) + REF_BASES_200.substring(0, 10);
        read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 1, readBases, readCigar);

        readContextCounter.processRead(read, 1, null);
        assertEquals(1, readContextCounter.readCounts().PartialCore);

        // and on the other side
        readCigar = "10M100N23M100N10M";

        readBases = REF_BASES_200.substring(0, 10) + altReadBases.substring(0, 23) + REF_BASES_200.substring(0, 10);
        read = createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 1, readBases, readCigar);

        readContextCounter.processRead(read, 1, null);
        assertEquals(2, readContextCounter.readCounts().PartialCore);
    }
}

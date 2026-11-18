package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.RegionTask.checkFragmentBoundaries;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.createSamRecord;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;
import static htsjdk.samtools.CigarOperator.S;

import java.util.List;

import com.hartwig.hmftools.common.bam.SupplementaryReadData;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class PositionEvidenceTest
{
    @Test
    public void equalsTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 100, "A", "C");
        PositionEvidence pe1 = new PositionEvidence("1", 100, "A", "C");
        PositionEvidence pe2 = new PositionEvidence("1", 100, "A", "G");
        PositionEvidence pe3 = new PositionEvidence("1", 100, "T", "C");
        PositionEvidence pe4 = new PositionEvidence("1", 101, "A", "C");
        PositionEvidence pe5 = new PositionEvidence("2", 100, "A", "C");
        assertEquals(pe0, pe1);
        Assert.assertNotEquals(pe0, pe2);
        Assert.assertNotEquals(pe0, pe3);
        Assert.assertNotEquals(pe0, pe4);
        Assert.assertNotEquals(pe0, pe5);

        PositionEvidence pe6 = new PositionEvidence("1", 100, "A", "C");
        assertEquals(pe0, pe6);
        pe6.ReadDepth = 10;
        assertEquals(pe0, pe6);
        pe6.ReadDepth = 1;
        assertEquals(pe0, pe6);
        pe6.RefSupport = 90;
        assertEquals(pe0, pe6);
        pe6.AltSupport = 10;
        assertEquals(pe0, pe6);
    }

    @Test
    public void hashCodeTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 100, "A", "C");
        PositionEvidence pe1 = new PositionEvidence("1", 100, "A", "C");
        assertEquals(pe0.hashCode(), pe1.hashCode());
        pe0.AltSupport = 10;
        pe0.RefSupport = 90;
        pe0.ReadDepth = 100;
        assertEquals(pe0.hashCode(), pe1.hashCode());
    }

    @Test
    public void convertToVafReadingTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 1000, "A", "C");
        pe0.AltSupport = 10;
        pe0.RefSupport = 90;
        pe0.ReadDepth = 100;
        VafReading converted = pe0.convertToVafReading();
        assertEquals("1", converted.chromosome());
        assertEquals(1000, converted.position());
        assertEquals(100, converted.readDepth());
        assertEquals(90, converted.refSupport());
        assertEquals(10, converted.altSupport());
    }

    @Test
    public void vafTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 1000, "A", "C");
        assertEquals(Double.NaN, pe0.vaf(), 0.0001);
        pe0.AltSupport = 10;
        pe0.RefSupport = 90;
        pe0.ReadDepth = 120;
        assertEquals(0.1, pe0.vaf(), 0.0001);
    }

    @Test
    public void symmetricVafTest()
    {
        PositionEvidence pe0 = new PositionEvidence("1", 1000, "A", "C");
        assertEquals(Double.NaN, pe0.symmetricVaf(), 0.0001);
        pe0.AltSupport = 90;
        pe0.RefSupport = 10;
        pe0.ReadDepth = 110;
        assertEquals(0.1, pe0.symmetricVaf(), 0.0001);

        PositionEvidence pe1 = new PositionEvidence("1", 1000, "A", "C");
        assertEquals(Double.NaN, pe1.symmetricVaf(), 0.0001);
        pe1.AltSupport = 20;
        pe1.RefSupport = 80;
        pe1.ReadDepth = 110;
        assertEquals(0.2, pe1.symmetricVaf(), 0.0001);
    }

    @Test
    public void testReadCigarState()
    {
        List<CigarElement> elements = List.of(new CigarElement(10, M));
        ReadCigarState state = new ReadCigarState(elements, 100, 10);
        assertEquals(0, state.CigarIndex);
        assertEquals(0, state.ElementIndex);
        assertEquals(0, state.ReadIndex);
        assertEquals(100, state.RefPosition);
        assertTrue(state.isValid());
        assertFalse(state.isDelete());
        assertFalse(state.isInsert());

        state.moveToPosition(109);
        assertEquals(0, state.CigarIndex);
        assertEquals(9, state.ElementIndex);
        assertEquals(9, state.ReadIndex);
        assertEquals(109, state.RefPosition);
        assertFalse(state.isDelete());
        assertFalse(state.isInsert());

        state.moveToPosition(110);
        assertFalse(state.isValid());

        // with delete and insert
        elements = List.of(
                new CigarElement(5, M),  // 100-104, index 0-4
                new CigarElement(2, D),  // 105-106 delete
                new CigarElement(5, M),  // 107-111, index 5-9
                new CigarElement(2, I),  // 111, index 10-11
                new CigarElement(5, M)); // 112-116, index 12-16

        state = new ReadCigarState(elements, 100, 20);
        assertEquals(0, state.CigarIndex); // pass the soft-clip
        assertEquals(0, state.ElementIndex);
        assertEquals(0, state.ReadIndex);
        assertEquals(100, state.RefPosition);

        state.moveToPosition(104);
        assertEquals(104, state.RefPosition);
        assertEquals(0, state.CigarIndex);
        assertEquals(4, state.ElementIndex);
        assertEquals(4, state.ReadIndex);
        assertFalse(state.isDelete());

        state.moveToPosition(105); // within the delete
        assertEquals(105, state.RefPosition);
        assertEquals(1, state.CigarIndex);
        assertEquals(0, state.ElementIndex);
        assertEquals(4, state.ReadIndex); // stays at last aligned base
        assertTrue(state.isDelete());

        state.moveToPosition(106); // within the delete
        assertEquals(106, state.RefPosition);
        assertEquals(1, state.CigarIndex);
        assertEquals(1, state.ElementIndex);
        assertEquals(4, state.ReadIndex);
        assertTrue(state.isDelete());

        state.moveToPosition(107);
        assertEquals(107, state.RefPosition);
        assertEquals(2, state.CigarIndex);
        assertEquals(0, state.ElementIndex);
        assertEquals(5, state.ReadIndex);
        assertFalse(state.isDelete());

        state.moveToPosition(111); // followed by insert
        assertEquals(111, state.RefPosition);
        assertEquals(2, state.CigarIndex);
        assertEquals(4, state.ElementIndex);
        assertEquals(9, state.ReadIndex);
        assertFalse(state.isDelete());
        assertTrue(state.isInsert());

        state.moveToPosition(112);
        assertEquals(112, state.RefPosition);
        assertEquals(4, state.CigarIndex);
        assertEquals(0, state.ElementIndex);
        assertEquals(12, state.ReadIndex);
        assertFalse(state.isInsert());

        state.moveToPosition(116);
        assertEquals(116, state.RefPosition);
        assertEquals(4, state.CigarIndex);
        assertEquals(4, state.ElementIndex);
        assertEquals(16, state.ReadIndex);

        // more complicated example with soft-clips
        elements = List.of(
                new CigarElement(5, S),  // index 0-4
                new CigarElement(10, M), // 100-109, index 5-14
                new CigarElement(10, I), // at 109, index 15-24
                new CigarElement(10, M), // 110-119, index 25-34
                new CigarElement(10, D), // 120-129, index stays at 34
                new CigarElement(10, M), // 130-139, index 35-44
                new CigarElement(5, S)); // index 45-49

        state = new ReadCigarState(elements, 100, 50);
        assertEquals(1, state.CigarIndex); // pass the soft-clip
        assertEquals(0, state.ElementIndex);
        assertEquals(5, state.ReadIndex);
        assertEquals(100, state.RefPosition);

        state.moveToPosition(105);
        assertEquals(10, state.ReadIndex);
        assertEquals(105, state.RefPosition);

        state.moveToPosition(109);
        assertEquals(14, state.ReadIndex);
        assertEquals(109, state.RefPosition);
        assertFalse(state.isDelete());
        assertTrue(state.isInsert());

        // first base past the insert
        state.moveToPosition(110);
        assertEquals(25, state.ReadIndex);
        assertEquals(110, state.RefPosition);
        assertFalse(state.isDelete());
        assertFalse(state.isInsert());

        state.moveToPosition(119);
        assertEquals(34, state.ReadIndex);
        assertEquals(119, state.RefPosition);
        assertFalse(state.isDelete());
        assertFalse(state.isInsert());

        state.moveToPosition(120);
        assertEquals(34, state.ReadIndex);
        assertEquals(120, state.RefPosition);
        assertTrue(state.isDelete());
        assertFalse(state.isInsert());

        state.moveToPosition(129);
        assertEquals(34, state.ReadIndex);
        assertEquals(129, state.RefPosition);
        assertTrue(state.isDelete());
        assertFalse(state.isInsert());

        state.moveToPosition(130);
        assertEquals(35, state.ReadIndex);
        assertEquals(130, state.RefPosition);
        assertFalse(state.isDelete());
        assertFalse(state.isInsert());

        state.moveToPosition(139);
        assertEquals(44, state.ReadIndex);
        assertEquals(139, state.RefPosition);

        state.moveToPosition(140);
        assertFalse(state.isValid());
    }

    @Test
    public void testFragmentReadSync()
    {
        String readId = "READ_001";
        String readBases = "A".repeat(20);

        SAMRecord read = createSamRecord(readId, CHR_1, 100, readBases, "20M", CHR_1, 110, false,
                false, null, true, "20M");

        int[] readCoords = checkFragmentBoundaries(read);
        assertEquals(100, readCoords[0]);
        assertEquals(119, readCoords[1]);

        // the mate only processes new / non-overlapping bases
        read = createSamRecord(readId, CHR_1, 110, readBases, "20M", CHR_1, 100, true,
                false, null, false, "20M");

        readCoords = checkFragmentBoundaries(read);
        assertEquals(120, readCoords[0]);
        assertEquals(129, readCoords[1]);

        // two reads overlapping at the start
        read = createSamRecord(readId, CHR_1, 100, readBases, "20M", CHR_1, 100, false,
                false, null, true, "25M");

        readCoords = checkFragmentBoundaries(read);
        assertEquals(100, readCoords[0]);
        assertEquals(119, readCoords[1]);

        // the mate still only gets the extras
        read = createSamRecord(readId, CHR_1, 100, readBases, "25M", CHR_1, 100, true,
                false, null, false, "20M");

        readCoords = checkFragmentBoundaries(read);
        assertEquals(120, readCoords[0]);
        assertEquals(124, readCoords[1]);

        // test 3' adapter trimming
        read = createSamRecord(readId, CHR_1, 100, readBases, "20M", CHR_1, 102, true,
                false, null, false, "20M");

        readCoords = checkFragmentBoundaries(read);
        assertEquals(102, readCoords[0]);
        assertEquals(119, readCoords[1]);

        read = createSamRecord(readId, CHR_1, 100, readBases, "20M", CHR_1, 110, true,
                false, null, false, "5S20M");

        readCoords = checkFragmentBoundaries(read);
        assertEquals(105, readCoords[0]);
        assertEquals(119, readCoords[1]);

        // on the other end
        read = createSamRecord(readId, CHR_1, 100, readBases, "20M", CHR_1, 90, false,
                false, null, true, "20M5S");

        readCoords = checkFragmentBoundaries(read);
        assertEquals(110, readCoords[0]); // due to overlapping bases
        assertEquals(114, readCoords[1]); // due to adapter trimming

        // finally conditions which skip any sync

        // BND
        read = createSamRecord(readId, CHR_1, 100, readBases, "20M", CHR_2, 110, false,
                false, null, true, "20M");

        assertNull(checkFragmentBoundaries(read));

        // INV
        read = createSamRecord(readId, CHR_1, 100, readBases, "20M", CHR_1, 110, false,
                false, null, false, "20M");

        assertNull(checkFragmentBoundaries(read));

        // unmapped mate
        read = createSamRecord(readId, CHR_1, 100, readBases, "20M", CHR_1, 110, false,
                false, null, true, "20M");
        read.setMateUnmappedFlag(true);

        assertNull(checkFragmentBoundaries(read));

        // beyond the other read
        read = createSamRecord(readId, CHR_1, 150, readBases, "20M", CHR_1, 110, false,
                false, null, true, "20M");

        assertNull(checkFragmentBoundaries(read));
    }
}

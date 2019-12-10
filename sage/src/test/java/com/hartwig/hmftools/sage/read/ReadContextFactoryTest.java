package com.hartwig.hmftools.sage.read;

import static org.junit.Assert.assertEquals;

import org.junit.Ignore;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

@Ignore
public class ReadContextFactoryTest {

    @Test
    public void testSimpleInsert() {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GAGGCTCATCTAGG";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("2M3I9M", readSequence);
        ReadContext victium = ReadContextFactory.createInsertContext("AGGC", 1000, 1, record, refBases);
        assertEquals("AGGCT", victium.centerBases());
    }

    @Test
    public void testInsertInRepeat() {
        String refSequence = "TGAAAAAAAATCT";
        String readSequence = "TGAAAAAAAAATCT";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("2M1I11M", readSequence);
        ReadContext victium = ReadContextFactory.createInsertContext("GA", 1000, 1, record, refBases);
        assertEquals("GAAAAAAAAAT", victium.centerBases());
    }

    @Test
    public void testInsertAtHomology() {
        String refSequence = "GATCATCTG";
        String readSequence = "GATCATCATCTG";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M3I8M", readSequence);
        ReadContext victium = ReadContextFactory.createInsertContext("ATCA", 1000, 1, record, refBases);
        assertEquals("GATCATCATCT", victium.centerBases());
    }

    @Test
    public void testInsertAtHomologyRepeat() {
        String refSequence = "GATCATCATCTG";
        String readSequence = "GATCATCATCATCTG";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M3I10M", readSequence);
        ReadContext victium = ReadContextFactory.createInsertContext("ATCA", 1000, 1, record, refBases);
        assertEquals("GATCATCATCATCT", victium.centerBases());
    }

    @Test
    public void testInsertAtHomologyWithAdditionalBases() {
        String refSequence = "ATGCGATCTTCC";
        String readSequence = "ATGCGATCAATCTTCC";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("5M4I7M", readSequence);
        ReadContext victium = ReadContextFactory.createInsertContext("GATCA", 1000, 4, record, refBases);
        assertEquals("GATCAA", victium.centerBases());
    }

    @Test
    public void testDeleteAtHomology() {
        String refSequence = "GATCGGATCGCTT";
        String readSequence = "GATCGCTT";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M5D7M", readSequence);
        ReadContext victium = ReadContextFactory.createDelContext("GATCGG", 1000, 0, record, refBases);
        assertEquals("GATCGC", victium.centerBases());
    }

    @Test
    public void testDeleteAtHomologyRepeat() {
        String refSequence = "GATCACCATCTG";
        String readSequence = "GATCATCTG";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M3D8M", readSequence);
        ReadContext victium = ReadContextFactory.createDelContext("ATCA", 1000, 1, record, refBases);
        assertEquals("GATCATCT", victium.centerBases());
    }

    @Test
    public void testDeleteAtRepeatInRef() {
        String refSequence = "GATCATCATCTG";
        String readSequence = "GATCATCTG";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M3D8M", readSequence);
        ReadContext victium = ReadContextFactory.createDelContext("ATCA", 1000, 1, record, refBases);
        assertEquals("GATCATCTG", victium.centerBases());
    }

    @Test
    public void testDeleteOneBase() {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GTCATCTAGG";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M1D9M", readSequence);
        ReadContext victium = ReadContextFactory.createDelContext("GA", 1000, 0, record, refBases);
        assertEquals("GT", victium.centerBases());
    }

    @Test
    public void testDeleteTwoBase() {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GCATCTAGG";
        IndexedBases refBases = new IndexedBases(1000, 0, refSequence.getBytes());

        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M2D8M", readSequence);
        ReadContext victium = ReadContextFactory.createDelContext("GAT", 1000, 0, record, refBases);
        assertEquals("GC", victium.centerBases());
    }

}

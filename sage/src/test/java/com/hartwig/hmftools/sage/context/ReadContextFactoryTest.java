package com.hartwig.hmftools.sage.context;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadContextFactoryTest {

    @Test
    public void testSimpleInsert() {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GAGGCTCATCTAGG";
        SAMRecord record = ReadContextDistanceTest.buildSamRecord("2M3I9M", readSequence);
        ReadContextImproved victium = ReadContextFactory.createInsertContext("AGGC", 1000, 1, record, 0, refSequence.getBytes());
        assertEquals("AGGCT", victium.centerBases());
    }

    @Test
    public void testInsertInRepeat() {
        String refSequence = "TGAAAAAAAATCT";
        String readSequence = "TGAAAAAAAAATCT";
        SAMRecord record = ReadContextDistanceTest.buildSamRecord("2M1I11M", readSequence);
        ReadContextImproved victium = ReadContextFactory.createInsertContext("GA", 1000, 1, record, 0, refSequence.getBytes());
        assertEquals("GAAAAAAAAAT", victium.centerBases());
    }

    @Test
    public void testInsertAtHomology() {
        String refSequence = "GATCATCTG";
        String readSequence = "GATCATCATCTG";
        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M3I8M", readSequence);
        ReadContextImproved victium = ReadContextFactory.createInsertContext("ATCA", 1000, 1, record, 0, refSequence.getBytes());
        assertEquals("GATCATCATCT", victium.centerBases());
    }

    @Test
    public void testInsertAtHomologyRepeat() {
        String refSequence = "GATCATCATCTG";
        String readSequence = "GATCATCATCATCTG";
        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M3I10M", readSequence);
        ReadContextImproved victium = ReadContextFactory.createInsertContext("ATCA", 1000, 1, record, 0, refSequence.getBytes());
        assertEquals("GATCATCATCATCT", victium.centerBases());
    }

    @Test
    public void testInsertAtHomologyWithAdditionalBases() {
        String refSequence =  "ATGCGATCTTCC";
        String readSequence = "ATGCGATCAATCTTCC";
        SAMRecord record = ReadContextDistanceTest.buildSamRecord("5M4I7M", readSequence);
        ReadContextImproved victium = ReadContextFactory.createInsertContext("GATCA", 1000, 4, record, 0, refSequence.getBytes());
        assertEquals("GATCAA", victium.centerBases());
    }

    @Test
    public void testDeleteAtHomology() {
        String refSequence =  "GATCGGATCGCTT";
        String readSequence = "GATCGCTT";
        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M5D7M", readSequence);
        ReadContextImproved victium = ReadContextFactory.createDelContext("GATCGG", 1000, 0, record, 1, refSequence.getBytes());
        assertEquals("GATCGC", victium.centerBases());
    }

    @Test
    public void testDeleteAtHomologyRepeat() {
        String refSequence = "GATCATCATCTG";
        String readSequence = "GATCATCTG";
        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M3D8M", readSequence);
        ReadContextImproved victium = ReadContextFactory.createDelContext("ATCA", 1000, 1, record, 1, refSequence.getBytes());
        assertEquals("GATCATCT", victium.centerBases());
    }

    @Test
    public void testDeleteOneBase() {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GTCATCTAGG";
        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M1D9M", readSequence);
        ReadContextImproved victium = ReadContextFactory.createDelContext("GA", 1000, 0, record, 0, refSequence.getBytes());
        assertEquals("GT", victium.centerBases());
    }

    @Test
    public void testDeleteTwoBase() {
        String refSequence = "GATCATCTAGG";
        String readSequence = "GCATCTAGG";
        SAMRecord record = ReadContextDistanceTest.buildSamRecord("1M2D8M", readSequence);
        ReadContextImproved victium = ReadContextFactory.createDelContext("GAT", 1000, 0, record, 0, refSequence.getBytes());
        assertEquals("GC", victium.centerBases());
    }

}

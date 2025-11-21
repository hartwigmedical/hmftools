package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.prep.JunctionsTest.REF_BASES;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.esvee.common.ReadIdTrimmer;

import org.checkerframework.dataflow.qual.Pure;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class MiscTest
{
    @Test
    public void testReadTrimmer()
    {
        ReadIdTrimmer trimmer = new ReadIdTrimmer(false);

        String readId = "READ_01";
        assertEquals(readId, trimmer.trim(readId));

        trimmer = new ReadIdTrimmer(true);

        // illumina style: A00260:251:HLYGFDSXY:1:1673:32280:4946
        readId = "A00260:251:HLYGFDSXY:1:1673:32280:4946";
        assertEquals("1:1673:32280:4946", trimmer.trim(readId));

        readId = "A00260:251:HLYGFDSXY:1:2:3:4";
        assertEquals("1:2:3:4", trimmer.trim(readId));

        // disable on a read with different delims or too short
        readId = "A0:25:HL:1:2:3:4";
        assertEquals(readId, trimmer.trim(readId));
        assertFalse(trimmer.enabled());

        trimmer = new ReadIdTrimmer(true); // reset

        readId = "A00260:251:HLYGFDSXY:1:1673:32280:4946";
        assertEquals("1:1673:32280:4946", trimmer.trim(readId));

        readId = "A00260:251:HLYGFDSX:1:2:3:4"; // diff delim position
        assertEquals(readId, trimmer.trim(readId));
        assertFalse(trimmer.enabled());

        // other formats are not currently supported
        trimmer = new ReadIdTrimmer(true); // reset

        readId = "011852_2-UGAv3-2-1333458495";
        assertEquals(readId, trimmer.trim(readId));
        assertFalse(trimmer.enabled());
    }

    @Test
    public void testDepthTracker()
    {
        DepthTracker depthTracker = new DepthTracker(new BaseRegion(1, 1000), 100);

        SAMRecord read = createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 1, REF_BASES.substring(0, 50), "50M");

        depthTracker.processRead(read);

        read = createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 10, REF_BASES.substring(0, 50), "50M");

        depthTracker.processRead(read);

        read = createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 120, REF_BASES.substring(0, 50), "50M");

        depthTracker.processRead(read);

        read = createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 180, REF_BASES.substring(0, 50), "50M");

        depthTracker.processRead(read);

        read = createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 1001, REF_BASES.substring(0, 50), "50M");

        depthTracker.processRead(read);

        read = createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 951, REF_BASES.substring(0, 50), "50M");

        depthTracker.processRead(read);

        assertEquals(1, depthTracker.calcDepth(1), 0.1);
        assertEquals(1, depthTracker.calcDepth(100), 0.1);
        assertEquals(0.7, depthTracker.calcDepth(101), 0.1);
        assertEquals(0.7, depthTracker.calcDepth(199), 0.1);
        assertEquals(0, depthTracker.calcDepth(301), 0.1);
        assertEquals(0.5, depthTracker.calcDepth(1000), 0.1);
        assertEquals(0, depthTracker.calcDepth(1001), 0.1);
    }
}

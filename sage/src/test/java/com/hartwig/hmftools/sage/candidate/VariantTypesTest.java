package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_MIN_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.RegionTaskTester;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.pipeline.RegionTask;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class VariantTypesTest
{
    private static final String TEST_REF_BASES = "X" // to cover the zero-position
            + "GCAGGAGAATCCCTTGAACCTGGGGGGAACCTGGGGGGGTGAGCTGAGAT"
            + "CATGCCATTGCACTCTAGCCTGGGCAACAAGAGTGAAACTCCGCCTCAAA"
            + "ACAAACAAACAAACAAACAAACAAACAAACAAACAAAAACCTCCAAAACA";
    // pos:    1         11        21        31        41       50 and repeats up to 150

    @Test
    public void testBasicVariants()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 200);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = TEST_REF_BASES;
        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500)); // need to cover the ref sequence buffer

        String readBases = refBases.substring(1, 21) + "A" + refBases.substring(22, 52);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 1, readBases, "50M");
        tester.TumorSamSlicer.ReadRecords.add(read1);

        readBases = refBases.substring(3, 21) + "A" + refBases.substring(22, 54);
        SAMRecord read2 = createSamRecord("READ_02", CHR_1, 3, readBases, "50M");
        tester.TumorSamSlicer.ReadRecords.add(read2);

        task.run();

        assertEquals(1, task.getVariants().size());
        SageVariant var = task.getVariants().get(0);
        assertEquals(21, var.position());
        assertEquals("T", var.ref());
        assertEquals("A", var.alt());
    }

    @Test
    public void testInsertWithRepeats()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 200);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = TEST_REF_BASES;
        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500)); // need to cover the ref sequence buffer

        String insertedBases = "GGGGGGGG";
        int readStartPos = 81;

        // extend the read to the end of the sequence so it's past the section of repeats
        String readBases = refBases.substring(readStartPos, 111) + insertedBases + refBases.substring(111, 151);

        // 71 = G, 67 = C, 65 = A, 84 = T
        SAMRecord read1 = createSamRecord("READ_01", CHR_1, readStartPos, readBases, "30M8I40M");
        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);

        task.run();

        assertEquals(1, task.getVariants().size());
        SageVariant var = task.getVariants().get(0);
        assertEquals(110, var.position());
        assertEquals("C", var.ref());
        assertEquals("C" + insertedBases, var.alt());
        assertFalse(var.candidate().readContext().hasIncompleteCore());
        assertEquals(40, var.candidate().readContext().indexedBases().coreLength());
    }

    @Test
    public void testInsertWithHomology()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 200);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = TEST_REF_BASES;
        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500)); // need to cover the ref sequence buffer

        String insertedBases = "AACCTGGGGGG";
        int readStartPos = 1;

        // extend the read to the end of the sequence so it's past the section of repeats
        String readBases = refBases.substring(readStartPos, 39) + insertedBases + refBases.substring(39, 81);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, readStartPos, readBases, "38M11I42M");
        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);

        task.run();

        assertEquals(1, task.getVariants().size());
        SageVariant var = task.getVariants().get(0);

        // no left alignment
        assertEquals(38, var.position());
        assertEquals("G", var.ref());
        assertEquals("G" + insertedBases, var.alt());
        assertFalse(var.candidate().readContext().hasIncompleteCore());
        assertEquals(39, var.candidate().readContext().indexedBases().coreLength());
    }

}

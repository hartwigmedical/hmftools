package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.SageConstants.SC_INSERT_MIN_FLANK_LENGTH;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;

import static junit.framework.TestCase.assertEquals;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.RegionTaskTester;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.pipeline.RegionTask;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class RegionTaskTest
{
    private static final String TEST_REF_BASES =
            "GCAGGAGAATCCCTTGAACCTGGGAGGCAGAGGTTACAGTGAGCTGAGAT"
          + "CATGCCATTGCACTCTAGCCTGGGCAACAAGAGTGAAACTCCGCCTCAAA"
          + "ACAAACAAACAAACAAACAAACAAACAAACAAACAAAAACCTCCAAAACA";
           //0         10        20        30        40       49

    @Test
    public void testBasicVariants()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 200);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = TEST_REF_BASES;
        tester.RefGenome.RefGenomeMap.put(CHR_1, "X" + refBases + generateRandomBases(1500)); // need to cover the ref sequence buffer

        String readBases = refBases.substring(0, 20) + "A" + refBases.substring(21, 51);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 1, readBases, "50M");
        tester.TumorSamSlicer.ReadRecords.add(read1);

        readBases = refBases.substring(2, 20) + "A" + refBases.substring(21, 53);
        SAMRecord read2 = createSamRecord("READ_02", CHR_1, 3, readBases, "50M");
        tester.TumorSamSlicer.ReadRecords.add(read2);

        task.run();

        assertEquals(1, task.getVariants().size());
        SageVariant var = task.getVariants().get(0);
        assertEquals(21, var.position());
        assertEquals("T", var.ref());
        assertEquals("A", var.alt());
    }

}

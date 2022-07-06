package com.hartwig.hmftools.sage.evidence;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.evidence.ReadContextCounter.RC_FULL;
import static com.hartwig.hmftools.sage.evidence.ReadContextCounter.RC_PARTIAL;
import static com.hartwig.hmftools.sage.evidence.ReadContextCounter.RC_REALIGNED;
import static com.hartwig.hmftools.sage.evidence.Realignment.realigned;
import static com.hartwig.hmftools.sage.evidence.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.evidence.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.evidence.RealignedType.NONE;
import static com.hartwig.hmftools.sage.evidence.RealignedType.SHORTENED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.RegionTaskTester;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.pipeline.RegionTask;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class RealignmentTest
{
    @Test
    public void testIndelRealignedLeft()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 450);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        // insert at position 288
        String refBases = "X" + generateRandomBases(187)
                + "CTACCCAAATCTTATAAAATGGCCCCACCCATATCTCCCTTCGCTGACTCTCTATTCGGACTCAGCCCGCCTGCACCCAGGTGAAATAAACAGCCTTGTT"
                + "GCACACACACACACACACACACACAAAGAATTCATTGCTCTTAGCTTTTGCACTGTTTCCTTTTTTAGACCTTCCTACATAAGTATTTGTGTATATGTCTGTATTT";

        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500));

        SAMRecord read0 = createSamRecord("READ_00", CHR_1, 257, // 24220, actual pos = 257
                "CCTGCACCCAGTTGAAATAAACAGCCTTGTTGCACACACACACACACACACACACACACACACACACACAAAGAATTCATTGCTCTTAGCTTTTGCACTGTTTCC"
                        + "TTTTTTAGACCTTCCTACATAAGTATTTGTGTATATGTCTGTATTT",
                "32M14I105M");

        tester.TumorSamSlicer.ReadRecords.add(read0);
        tester.TumorSamSlicer.ReadRecords.add(read0);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 188, // 24220, actual pos =
                "CTACCCAAATCTTATAAAATGGCCCCACCCATATCTCCCTTCGCTGACTCTCTATTCGGACTCAGCCCGCCTGCACCCAGGTGAAATAAACAGCCTTGTTGC"
                        + "ACACACACACACACACACACACACACACACACACACAAAGAATTCATTG",
                "125M26S");

        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);

        task.run();

        SageVariant var1 = task.getVariants().stream().filter(x -> x.position() == 288).findFirst().orElse(null);
        SageVariant var2 = task.getVariants().stream().filter(x -> x.position() == 312).findFirst().orElse(null);
        assertNotNull(var1);
        assertNotNull(var2);
        assertEquals(4, var1.tumorReadCounters().get(0).counts()[RC_FULL]);
        assertEquals(2, var2.tumorReadCounters().get(0).counts()[RC_FULL]);
        assertEquals(2, var2.tumorReadCounters().get(0).counts()[RC_REALIGNED]);
    }

    @Test
    public void testIndelRealignedRight()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 450);

        RegionTaskTester tester = new RegionTaskTester();

        RegionTask task = tester.createRegionTask(region);

        String refBases = "X"
                + "AATTTGAACCTAATTTTTTTTTTTCTTCTGCACTAACATGCCTGTTGAACCATTTGGACTTAACTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCATTCTTTTGCATGTAGATAT"
                + "CCTTTTCCCAGCACCATTCGTTGAATGGAGACTATTCTTTCCCCACTGAATAGTCTTGGTACCCTCTTTGAAAATCAATTGATGATAAATAGATGTGTTTATTTCTGAACTCTCCATTT"
                + "TATTCCATTGACCTATATCTCTCCTTATGCCAGTTTTTATTACTGTGCAGTTTTGATTACTAC";

        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1500));

        SAMRecord read0 = createSamRecord("READ_00", CHR_1, 28, // 11490
                "CTGCACTAACATGCCTGTTGAACCATTTGGACTTAACTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCATT"
                        + "CTTTTGCATGTAGATATCCTTTTCCCAGCACCATTCGTTGAAT",
                "36M34I81M");

        tester.TumorSamSlicer.ReadRecords.add(read0);
        tester.TumorSamSlicer.ReadRecords.add(read0);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 14, // 15718
                "ATTTTTTTTTTCTTCTGCACTAACATGCCTGTTGAACCATTTGGACTTAACTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCTTTTGTGCATGGTGTGAAATAGGTG"
                        + "CCCAGCCTCATTCTTTTGCATGTAGATATCCTTTTCCCAGC",
                "85M66S");

        tester.TumorSamSlicer.ReadRecords.add(read1);
        tester.TumorSamSlicer.ReadRecords.add(read1);

        SAMRecord read2 = createSamRecord("READ_02", CHR_1, 64, // 22780
                "TTGGACTTAACTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCTTTTGTGCATGGTGTGAAATAGGTGCCCAGCCTCATTCTTTTGCATGTAGATATCCTTTTCCCAG"
                        + "CACCATTCGTTGAATGGAGACTATTCTTTCCCCACTGAATA",
                "44S107M");

        tester.TumorSamSlicer.ReadRecords.add(read2);
        tester.TumorSamSlicer.ReadRecords.add(read2);

        task.run();

        SageVariant var1 = task.getVariants().stream().filter(x -> x.position() == 63).findFirst().orElse(null);
        SageVariant var2 = task.getVariants().stream().filter(x -> x.position() == 98).findFirst().orElse(null);
        assertNotNull(var1);
        assertNotNull(var2);
        assertEquals(4, var1.tumorReadCounters().get(0).counts()[RC_FULL]);
        assertEquals(2, var1.tumorReadCounters().get(0).counts()[RC_PARTIAL]);

        assertEquals(2, var2.tumorReadCounters().get(0).counts()[RC_FULL]);
        // assertEquals(2, var1.tumorReadCounters().get(0).counts()[RC_PARTIAL]);
        // assertEquals(4, var2.tumorReadCounters().get(0).counts()[RC_REALIGNED]);
    }

    @Test
    public void testRealignedTooShort()
    {
        String sequence = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTT";
        String truncatedAtEnd = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTT";
        String truncatedAtStart = "GAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTT";
        int startIndex = 0;
        int endIndex = sequence.length() - 1;

        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex, sequence.getBytes(), 0));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex, truncatedAtEnd.getBytes(), 0));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), -2, truncatedAtStart.getBytes(), 0));
    }

    @Test
    public void testRealigned()
    {
        String sequence = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTTATTATTG";
        int startIndex = 3;
        int endIndex = startIndex + 55;

        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex, sequence.getBytes(), 0));

        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 1, sequence.getBytes(), 0));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 1, sequence.getBytes(), 0));

        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 2, sequence.getBytes(), 1));
        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 1, sequence.getBytes(), 1));
        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 1, sequence.getBytes(), 1));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 2, sequence.getBytes(), 1));

        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 2, sequence.getBytes(), 2));
        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 2, sequence.getBytes(), 2));
    }

    @Test
    public void testPolyA()
    {
        String shorter = "GATCAAAAAAAAAGATC";
        String ref = "GATCAAAAAAAAAAGATC";
        String longer = "GATCAAAAAAAAAAAGATC";

        int startIndex = 0;
        int endIndex = ref.length() - 1;

        assertRealigned(EXACT, realigned(startIndex, endIndex, ref.getBytes(), startIndex, ref.getBytes(), 10));
        assertRealigned(SHORTENED, 9, realigned(startIndex, endIndex, ref.getBytes(), startIndex, shorter.getBytes(), 10));
        assertRealigned(LENGTHENED, 11, realigned(startIndex, endIndex, ref.getBytes(), startIndex, longer.getBytes(), 10));
    }

    @Test
    public void testDiNucleotideRepeat()
    {
        String shorter = "GATCATATATATGATC";
        String ref = "GATCATATATATATGATC";
        String longer = "GATCATATATATATATGATC";

        int startIndex = 0;
        int endIndex = ref.length() - 1;

        assertRealigned(EXACT, realigned(startIndex, endIndex, ref.getBytes(), startIndex, ref.getBytes(), 10));
        assertRealigned(SHORTENED, 4, realigned(startIndex, endIndex, ref.getBytes(), startIndex, shorter.getBytes(), 10));
        assertRealigned(LENGTHENED, 6, realigned(startIndex, endIndex, ref.getBytes(), startIndex, longer.getBytes(), 10));
    }

    private static void assertRealigned(RealignedType expectedType, int expectedCount, RealignedContext context)
    {
        assertEquals(expectedCount, context.RepeatCount);
        assertEquals(expectedType, context.Type);
    }

    private static void assertRealigned(RealignedType expectedType, RealignedContext context)
    {
        assertEquals(expectedType, context.Type);
    }
}

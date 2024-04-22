package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.common.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.RegionTaskTester;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.VariantReadContext;
import com.hartwig.hmftools.sage.pipeline.RegionTask;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class CandidateCreationTest
{
    @Test
    public void testSoftClipInsert()
    {
        String refBaseStr = generateRandomBases(100);

        RefSequence refBases = new RefSequence(100, refBaseStr.getBytes());

        String insertBases = "AAAAA";

        // first an insert on the right
        int scStartRefIndex = 51;
        String scBases = insertBases + refBaseStr.substring(scStartRefIndex, 71);
        String readBases = refBaseStr.substring(20, scStartRefIndex) + scBases;
        int scReadIndex = 31;

        AltRead altRead = RefContextConsumer.processSoftClip(
                120, 150, readBases, scBases.length(), scReadIndex, refBases, false);

        assertNotNull(altRead);
        assertEquals("G", altRead.Ref);
        assertEquals("G" + insertBases, altRead.Alt);

        // then on the left
        scStartRefIndex = 0;
        scBases = refBaseStr.substring(scStartRefIndex, 20) + insertBases;
        readBases = scBases + refBaseStr.substring(20, 51);
        scReadIndex = 0;

        altRead = RefContextConsumer.processSoftClip(
                120, 150, readBases, scBases.length(), scReadIndex, refBases, true);

        assertNotNull(altRead);
        assertEquals("T", altRead.Ref);
        assertEquals("T" + insertBases, altRead.Alt);
    }

    @Test
    public void testSoftClipInsertProdExample()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 150);

        RegionTaskTester tester = new RegionTaskTester();

        tester.PanelRegions.add(new BaseRegion(1, 1000));

        RegionTask task = tester.createRegionTask(region);

        // bases from 21,974,750 -> 950
        String refBases = "TCTACCCGACCCCGGGCCGCGGCCGTGGCCAGCCAGTCAGCCGAAGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCCGCTGCCTGCTCTCCCC"
                + "CTCTCCGCAGCCGCCGAGCGCACGCGGTCCGCCCCACCCTCTGGTGACCAGCCAGCCCCTCCTCTTTCTTCCTCCGGTGCTGGCGGAAGAG";

        tester.RefGenome.RefGenomeMap.put(CHR_1, refBases + generateRandomBases(1100)); // need to cover the ref sequence buffer

        List<SAMRecord> reads = Lists.newArrayList(
                createSamRecord(
                        READ_ID_GENERATOR.nextId(), CHR_1,  45,
                        "TGGCCAGCCAGTCAGCCGAAGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCCGCTGCCTGAAA",
                        "44S57M"),
                createSamRecord(
                        READ_ID_GENERATOR.nextId(), CHR_1,  45,
                        "GCCAGCCAGTCAGCCGAAGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCCGCACCCTGCTAAA",
                        "42S59M"),
                createSamRecord(
                        READ_ID_GENERATOR.nextId(), CHR_1,  45,
                        "CGCCAGGCAGCCGAAGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCCGCTGCCTGCTCTCAAA",
                        "39S62M"),
                createSamRecord(
                        READ_ID_GENERATOR.nextId(), CHR_1,  45,
                        "AGCCGAAGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCCGCTGCCTGCTCTCCCCCTCTCAAA",
                        "31S70M"));

        reads.get(1).setFirstOfPairFlag(false);
        reads.get(3).setFirstOfPairFlag(false);

        tester.TumorSamSlicer.ReadRecords.addAll(reads);
        tester.TumorSamSlicer.ReadRecords.addAll(reads); // repeat to get over qual thresholds
        tester.TumorSamSlicer.ReadRecords.addAll(reads);

        task.run();

        SageVariant var = task.getVariants().stream().filter(x -> x.position() == 44 && x.isIndel()).findFirst().orElse(null);

        Assert.assertNotNull(var);

        assertEquals(44, var.position());
        assertEquals("A", var.ref());
        assertEquals("AGGCTCCATGCTGCTCCCCGCCGCC", var.alt());

        VariantReadContext readContext = var.readContext();
        assertEquals(10, readContext.CoreIndexStart);
        assertEquals(11, readContext.VarReadIndex);
        assertEquals(11, readContext.VarReadIndex);
        assertEquals("36S57M", readContext.readCigar());
        assertEquals(48, readContext.Homology.Length);
    }
}

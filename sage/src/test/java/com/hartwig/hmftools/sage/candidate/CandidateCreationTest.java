package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;
import static com.hartwig.hmftools.sage.vcf.VariantVCF.DEDUP_INDEL_FILTER;

import static org.junit.Assert.assertTrue;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.RegionTaskTester;
import com.hartwig.hmftools.sage.common.SageVariant;
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

        IndexedBases refBases = new IndexedBases(100, 0, refBaseStr.getBytes());

        String insertBases = "AAAAA";

        // first an insert on the right
        int scStartRefIndex = 51;
        String scBases = insertBases + refBaseStr.substring(scStartRefIndex, 71);
        String readBases = refBaseStr.substring(20, scStartRefIndex) + scBases;
        int scReadIndex = 31;

        AltRead altRead = RefContextConsumer.processSoftClip(
                120, 150, readBases, scBases.length(), scReadIndex, refBases, false);

        assertNotNull(altRead);
        assertEquals("T", altRead.Ref);
        assertEquals("T" + insertBases, altRead.Alt);

        // then on the left
        scStartRefIndex = 0;
        scBases = refBaseStr.substring(scStartRefIndex, 20) + insertBases;
        readBases = scBases + refBaseStr.substring(20, 51);
        scReadIndex = 0;

        altRead = RefContextConsumer.processSoftClip(
                120, 150, readBases, scBases.length(), scReadIndex, refBases, true);

        assertNotNull(altRead);
        assertEquals("C", altRead.Ref);
        assertEquals("C" + insertBases, altRead.Alt);
    }

    @Test
    public void testSoftClipInsertProd1()
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
                        "36996", CHR_1,  45,
                        "TGGCCAGCCAGTCAGCCGAAGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCCGCTGCCTG",
                        "44S57M"),
                createSamRecord(
                        "8077", CHR_1,  45,
                        "GCCAGCCAGTCAGCCGAAGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCCGCACCCTGCT",
                        "42S59M"),
                createSamRecord(
                        "19351", CHR_1,  45,
                        "CGCCAGGCAGCCGAAGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCCGCTGCCTGCTCTC",
                        "39S62M"),
                createSamRecord(
                        "24001", CHR_1,  45,
                        "AGCCGAAGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCGGCTCCATGCTGCTCCCCGCCGCCCGCTGCCTGCTCTCCCCCTCTC",
                        "31S70M"));

        reads.get(1).setFirstOfPairFlag(false);
        reads.get(3).setFirstOfPairFlag(false);

        tester.TumorSamSlicer.ReadRecords.addAll(reads);
        tester.TumorSamSlicer.ReadRecords.addAll(reads); // repeat to get over qual thresholds
        tester.TumorSamSlicer.ReadRecords.addAll(reads);

        task.run();

        SageVariant var = task.getVariants().stream().filter(x -> x.position() == 44 && x.isIndel()).findFirst().orElse(null);

        Assert.assertNotNull(var);
        assertEquals(12, var.tumorReadCounters().get(0).softClipInsertSupport());
    }
}

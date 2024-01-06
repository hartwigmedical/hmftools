package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_MISMATCH_READS;
import static com.hartwig.hmftools.esvee.SvConstants.PRIMARY_ASSEMBLY_MIN_MISMATCH_TOTAL_QUAL;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.common.AssemblyUtils.purgeLowSupport;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.esvee.common.AssemblyMismatchSplitter;
import com.hartwig.hmftools.esvee.common.AssemblySequence;
import com.hartwig.hmftools.esvee.common.BaseMismatch;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

import org.junit.Test;

public class AssemblySequenceTest
{
    @Test
    public void testBuildJunctionSequence()
    {
        Read read1 = createSamRecord("READ_01", 10, REF_BASES.substring(10, 30) + "AACCGG", "20M6S");
        Read read2 = createSamRecord("READ_02", 9, REF_BASES.substring(9, 30)   + "ACCCG", "21M5S");
        Read read3 = createSamRecord("READ_03", 15, REF_BASES.substring(15, 30) + "AATCGGTT", "15M8S");
        Read read4 = createSamRecord("READ_04", 10, REF_BASES.substring(10, 30) + "AATCGGG", "20M7S");

        Junction junction = new Junction(CHR_1, 29, POS_STRAND);

        AssemblySequence junctionSequence = new AssemblySequence(junction, read1, 29, 37);
        junctionSequence.addRead(read2, true);
        junctionSequence.addRead(read3, true);
        junctionSequence.addRead(read4, true);

        List<BaseMismatch> baseMismatches = junctionSequence.mismatches().allBaseMismatches();
        assertEquals(3, baseMismatches.size());

        assertTrue(purgeLowSupport(junctionSequence, PRIMARY_ASSEMBLY_MIN_MISMATCH_READS, PRIMARY_ASSEMBLY_MIN_MISMATCH_TOTAL_QUAL));

        AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(junctionSequence);
        List<AssemblySequence> allSequences = splitter.splitOnMismatches(5);
        assertEquals(2, allSequences.size());

        assertEquals(1, allSequences.get(0).supportCount());
        assertEquals(2, allSequences.get(1).supportCount());
    }
}

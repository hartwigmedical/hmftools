package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES;
import static com.hartwig.hmftools.esvee.TestUtils.cloneRead;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.esvee.common.AssemblyMismatchSplitter;
import com.hartwig.hmftools.esvee.common.JunctionAssembly;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

import org.junit.Test;

public class JunctionAssemblyTest
{
    @Test
    public void testBuildJunctionSequence()
    {
        String refBases = REF_BASES.substring(0, 20) + "TTGGCCAATT";

        //   pos           10        20         30
        //   pos 012345678901234567890123456789 0123456789
        //       GATCGATCGATCGATCGATCTTGGCCAATT AACCGGGG (first sequence)
        // 1st-asm index  012345678901234567890 123456

        //       GATCGATCGATCGATCGATCTTGGCCAATT AATCGGTT (2nd sequence)
        // 2nd-asm index  012345678901234567890 12345678

        Junction junction = new Junction(CHR_1, 29, POS_STRAND);

        // read 1 defines the first sequence
        Read read1 = createSamRecord("READ_01", 10, refBases.substring(10, 30) + "AACCGGGG", "20M8S");
        Read read1b = cloneRead(read1, read1.getName() + "b");

        // read 2 has one mismatch on the first sequence and has a ref base mismatch
        String readRefBases = REF_BASES.substring(0, 20) + "TGGGCCAATT";
        Read read2 = createSamRecord("READ_02", 9, readRefBases.substring(9, 30) + "ACCCGG", "21M6S");

        // read 3 defines the second sequence
        Read read3 = createSamRecord("READ_03", 15, refBases.substring(15, 30) + "AATCGGTT", "15M8S");

        // read 4 matches has mismatches against both sequences, but only one with the second
        Read read4 = createSamRecord("READ_04", 10, refBases.substring(10, 30) + "AATCGGA", "20M7S");

        // read 5 matches the second sequence and also has a ref base mismatch
        readRefBases = REF_BASES.substring(0, 20) + "TTGCCCAATT";
        Read read5 = createSamRecord("READ_05", 10, readRefBases.substring(10, 30) + "AATCGGT", "20M7S");

        // has 1 mismatch against each so will support both
        Read read6 = createSamRecord("READ_06", 10, refBases.substring(10, 30) + "AATCGGG", "20M7S");

        JunctionAssembly junctionSequence = new JunctionAssembly(junction, read1, 29, 37);
        junctionSequence.addRead(read2, true);
        junctionSequence.addRead(read1b, true);
        junctionSequence.addRead(read3, true);
        junctionSequence.addRead(read4, true);
        junctionSequence.addRead(read5, true);
        junctionSequence.addRead(read6, true);

        assertEquals(5, junctionSequence.mismatches().allBaseMismatches().size());
        assertEquals(7, junctionSequence.supportCount());
        assertEquals(4, junctionSequence.mismatches().positionCount());

        AssemblyMismatchSplitter splitter = new AssemblyMismatchSplitter(junctionSequence);
        List<JunctionAssembly> allSequences = splitter.splitOnMismatches(5);
        assertEquals(2, allSequences.size());

        assertEquals(4, allSequences.get(0).supportCount());
        assertEquals(4, allSequences.get(1).supportCount());

        allSequences.forEach(x -> x.expandReferenceBases());

        // for now support count is duplicated across reads for ref bases and post-junction bases
        assertEquals(4, allSequences.get(0).supportCount());
        assertEquals(0, allSequences.get(0).mismatches().positionCount());
        assertEquals(4, allSequences.get(1).supportCount());
        assertEquals(0, allSequences.get(1).mismatches().positionCount());

    }
}

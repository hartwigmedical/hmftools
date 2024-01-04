package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Strand.NEG_STRAND;
import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES;
import static com.hartwig.hmftools.esvee.TestUtils.TEST_READ_ID;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.assembly.PrimaryAssembler.realignForJunction;
import static com.hartwig.hmftools.esvee.read.Read.INVALID_INDEX;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.esvee.common.AssemblySequence;
import com.hartwig.hmftools.esvee.common.BaseMismatch;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

import org.junit.Test;

public class PrimaryAssemblyTest
{
    @Test
    public void testReadIndelToSoftClips()
    {
        Junction fwdJunction = new Junction(CHR_1, 30, POS_STRAND);

        Read read = createSamRecord(TEST_READ_ID, 20, REF_BASES.substring(20, 40), "10M2I8M");
        assertEquals(12, read.getReadIndexAtReferencePosition(fwdJunction.Position)); // being the first base after the insert

        Read realignedRead = realignForJunction(read, fwdJunction);
        assertEquals(10, realignedRead.getReadIndexAtReferencePosition(fwdJunction.Position, true));

        assertEquals(29, realignedRead.getAlignmentEnd());
        assertEquals(39, realignedRead.getUnclippedEnd());

        Junction revJunction = new Junction(CHR_1, 30, NEG_STRAND);

        read = createSamRecord(
                TEST_READ_ID, 23, REF_BASES.substring(23, 31) + "AA" + REF_BASES.substring(31, 40), "8M2I10M");

        assertEquals(7, read.getReadIndexAtReferencePosition(fwdJunction.Position, true));

        realignedRead = realignForJunction(read, revJunction);
        assertEquals(9, realignedRead.getReadIndexAtReferencePosition(fwdJunction.Position, true));

        assertEquals(31, realignedRead.getAlignmentStart());
        assertEquals(31, realignedRead.getAlignmentStart());
        assertEquals(21, realignedRead.getUnclippedStart());
    }

    @Test
    public void testInitialAssemblySequence()
    {
        Read read1 = createSamRecord("READ_01", 10, REF_BASES.substring(10, 30) + "AACCGG", "20M6S");
        Read read2 = createSamRecord("READ_02", 9, REF_BASES.substring(9, 30)   + "ACCCG", "21M5S");
        Read read3 = createSamRecord("READ_03", 15, REF_BASES.substring(15, 30) + "AATCGGTT", "15M8S");
        Read read4 = createSamRecord("READ_04", 10, REF_BASES.substring(10, 30) + "AATCGGG", "20M7S");

        Junction junction = new Junction(CHR_1, 29, POS_STRAND);

        AssemblySequence assemblySequence = new AssemblySequence(junction, read1, 29, 37);
        assemblySequence.tryAddRead(read2);
        assemblySequence.tryAddRead(read3);
        assemblySequence.tryAddRead(read4);

        List<BaseMismatch> baseMismatches = assemblySequence.baseMismatches();
        assertEquals(3, baseMismatches.size());
    }

}

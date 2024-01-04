package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Strand.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.esvee.TestUtils.createSAMRecord;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.esvee.common.AssemblySequence;
import com.hartwig.hmftools.esvee.common.BaseMismatch;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

import org.junit.Test;

public class AssemblySequenceTest
{
    private static final String REF_BASES = generateRandomBases(100);

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

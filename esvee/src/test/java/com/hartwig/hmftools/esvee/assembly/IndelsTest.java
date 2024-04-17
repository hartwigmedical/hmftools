package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.read.Read;

import org.junit.Test;

public class IndelsTest
{
    @Test
    public void testLongDeleteAssemblies()
    {
        Junction posJunction = new Junction(CHR_1, 50, POS_ORIENT, false, true, false);

        // first a basic assembly with all reads agreeing
        String readBases = REF_BASES_200.substring(11, 51) + REF_BASES_200.substring(100, 140);
        Read read1 = createRead(READ_ID_GENERATOR.nextId(), 11, readBases, "40M49D40M");

        readBases = REF_BASES_200.substring(21, 51) + REF_BASES_200.substring(100, 150);
        Read read2 = createRead(READ_ID_GENERATOR.nextId(), 21, readBases, "30M49D50M");

        // other reads will soft-clip at the junctions
        readBases = REF_BASES_200.substring(11, 51) + REF_BASES_200.substring(100, 120);
        Read read3 = createRead(READ_ID_GENERATOR.nextId(), 11, readBases, "40M20S");

        List<Read> reads = List.of(read1, read2, read3);

        JunctionAssembler junctionAssembler = new JunctionAssembler(posJunction);
        List<JunctionAssembly> assemblies = junctionAssembler.processJunction(reads);
        assertEquals(1, assemblies.size());
        JunctionAssembly assembly = assemblies.get(0);
        assertEquals(3, assembly.supportCount());
        assertEquals(0, assembly.mismatchReadCount());

        // test the other side
        Junction negJunction = new Junction(CHR_1, 100, NEG_ORIENT, false, true, false);

        readBases = REF_BASES_200.substring(31, 51) + REF_BASES_200.substring(100, 160);
        read3 = createRead(READ_ID_GENERATOR.nextId(), 100, readBases, "20S60M");

        reads = List.of(read1, read2, read3);

        junctionAssembler = new JunctionAssembler(negJunction);
        assemblies = junctionAssembler.processJunction(reads);
        assertEquals(1, assemblies.size());
        assembly = assemblies.get(0);
        assertEquals(3, assembly.supportCount());
    }
}

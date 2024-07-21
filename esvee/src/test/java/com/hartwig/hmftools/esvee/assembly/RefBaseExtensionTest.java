package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.REF_GENOME;
import static com.hartwig.hmftools.esvee.TestUtils.createRead;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.extendRefBases;
import static com.hartwig.hmftools.esvee.assembly.RefBaseExtender.trimAssemblyRefBases;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.read.Read;

import org.junit.Test;

public class RefBaseExtensionTest
{
    @Test
    public void testBasicRefBaseExtension()
    {
        REF_GENOME.RefGenomeMap.put(CHR_1, REF_BASES_200);

        Junction posJunction = new Junction(CHR_1, 150, FORWARD);

        String initialRefBases = REF_BASES_200.substring(130, 151);
        String assemblyBases = initialRefBases + REF_BASES_200.substring(180);
        byte[] baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());
        int junctionIndex = initialRefBases.length() - 1;

        JunctionAssembly assembly = new JunctionAssembly(posJunction, assemblyBases.getBytes(), baseQuals, junctionIndex);
        assertEquals(130, assembly.refBasePosition());

        // create a collection of junction mate and discordant reads to test as support
        Read mate1 = createRead(READ_ID_GENERATOR.nextId(), 121, REF_BASES_200.substring(121, 151), "30M");
        mate1.markJunctionMate();

        String mate2Bases = REF_BASES_200.substring(121, 140) + "A" + REF_BASES_200.substring(141, 151);
        Read mate2 = createRead(READ_ID_GENERATOR.nextId(), 121, mate2Bases, "30M");
        mate2.markJunctionMate();

        // too many mismatches
        String disc1Bases = REF_BASES_200.substring(100, 115) + "TTT" + REF_BASES_200.substring(118, 130);
        Read disc1 = createRead(READ_ID_GENERATOR.nextId(), 100, disc1Bases, "30M");

        // supports the ref bases
        Read disc2 = createRead(READ_ID_GENERATOR.nextId(), 90, REF_BASES_200.substring(90, 120), "30M");

        // supports the ref bases but too large a gap so will be cleaned up later
        Read disc3 = createRead(READ_ID_GENERATOR.nextId(), 30, REF_BASES_200.substring(30, 60), "30M");

        List<Read> candidateSupport = Lists.newArrayList(mate1, mate2, disc1, disc2, disc3);

        extendRefBases(assembly, candidateSupport, REF_GENOME, false);

        assertEquals(30, assembly.refBasePosition());
        assertEquals(4, assembly.supportCount());

        trimAssemblyRefBases(assembly, 20);

        assertEquals(90, assembly.refBasePosition());
        assertEquals(REF_BASES_200.substring(90, 151), assembly.formRefBaseSequence());


        // repeat for -ve junction
        Junction negJunction = new Junction(CHR_1, 50, REVERSE);

        initialRefBases = REF_BASES_200.substring(50, 71);
        assemblyBases = REF_BASES_200.substring(0, 20) + initialRefBases;
        baseQuals = SamRecordTestUtils.buildDefaultBaseQuals(assemblyBases.length());
        junctionIndex = 20;

        assembly = new JunctionAssembly(negJunction, assemblyBases.getBytes(), baseQuals, junctionIndex);
        assertEquals(70, assembly.refBasePosition());

        // create a collection of junction mate and discordant reads to test as support
        mate1 = createRead(READ_ID_GENERATOR.nextId(), 50, REF_BASES_200.substring(50, 80), "30M");
        mate1.markJunctionMate();

        mate2Bases = REF_BASES_200.substring(60, 75) + "A" + REF_BASES_200.substring(76, 90);
        mate2 = createRead(READ_ID_GENERATOR.nextId(), 60, mate2Bases, "30M");
        mate2.markJunctionMate();

        // supports the ref bases but too large a gap so will be cleaned up later
        disc1 = createRead(READ_ID_GENERATOR.nextId(), 120, REF_BASES_200.substring(120, 150), "30M");

        // too many mismatches
        String disc2Bases = REF_BASES_200.substring(100, 115) + "TTT" + REF_BASES_200.substring(118, 160);
        disc2 = createRead(READ_ID_GENERATOR.nextId(), 100, disc2Bases, "60M");

        candidateSupport = Lists.newArrayList(mate1, mate2, disc1, disc2);

        extendRefBases(assembly, candidateSupport, REF_GENOME, false);

        assertEquals(159, assembly.refBasePosition());
        assertEquals(3, assembly.supportCount());

        trimAssemblyRefBases(assembly, 20);

        assertEquals(89, assembly.refBasePosition());
        assertEquals(REF_BASES_200.substring(50, 90), assembly.formRefBaseSequence());
    }
}

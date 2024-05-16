package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.TestUtils.createAssembly;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.esvee.assembly.phase.JunctionExtender;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;

import org.junit.Test;

public class JunctionExtenderTest
{
    private static final String EXT_REF_BASES = REF_BASES_200 + REF_BASES_200;

    @Test
    public void testJunctionExtensions()
    {
        String junctionBases = REF_BASES_200.substring(0, 150);
        JunctionAssembly assembly = createAssembly(CHR_1, 200, FORWARD, junctionBases, 99);

        String otherAssemblyBases = REF_BASES_200.substring(120, 180);
        JunctionAssembly otherAssembly = createAssembly(CHR_2, 200, FORWARD, otherAssemblyBases, 20);

        JunctionExtender junctionExtender = new JunctionExtender(assembly, List.of(otherAssembly), Collections.emptyList());
        junctionExtender.extendAssembly();



    }


}

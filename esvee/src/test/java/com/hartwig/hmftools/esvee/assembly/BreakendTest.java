package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_200;
import static com.hartwig.hmftools.esvee.assembly.AssemblyTestUtils.createAssemblyAlignment;

import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.esvee.alignment.AssemblyAlignment;
import com.hartwig.hmftools.esvee.alignment.Breakend;
import com.hartwig.hmftools.esvee.alignment.Deduplication;
import com.hartwig.hmftools.esvee.alignment.HomologyData;
import com.hartwig.hmftools.esvee.common.FilterType;

import org.junit.Test;

public class BreakendTest
{
    @Test
    public void testDeduplication()
    {
        MockRefGenome refGenome = new MockRefGenome();
        refGenome.RefGenomeMap.put(CHR_1, REF_BASES_200);
        refGenome.RefGenomeMap.put(CHR_2, REF_BASES_200);

        AssemblyAlignment assemblyAlignment = createAssemblyAlignment(
                refGenome, CHR_1, 100, FORWARD, CHR_2, 50, REVERSE, "", "");

        Breakend breakend1 = new Breakend(assemblyAlignment, CHR_1, 100, FORWARD, "", null);
        Breakend breakend2 = new Breakend(assemblyAlignment, CHR_1, 100, REVERSE, "", null);

        Breakend breakend3 = new Breakend(assemblyAlignment, CHR_1, 100, FORWARD, "", null);

        HomologyData homology = new HomologyData("AAA", -3, 3, -3, 3);
        Breakend breakend4 = new Breakend(assemblyAlignment, CHR_1, 102, FORWARD, "", homology);
        Breakend breakend5 = new Breakend(assemblyAlignment, CHR_1, 104, FORWARD, "", homology);

        List<Breakend> breakends = Lists.newArrayList(breakend1, breakend2, breakend3, breakend4, breakend5);
        Collections.sort(breakends);

        Deduplication.deduplicateBreakends(breakends);

        assertTrue(breakend1.passing());
        assertTrue(breakend2.passing());
        assertTrue(breakend3.filters().contains(FilterType.DUPLICATE));
        assertTrue(breakend4.filters().contains(FilterType.DUPLICATE));
        assertTrue(breakend5.passing());
    }
}

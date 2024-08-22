package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.caller.CallerTestUtils.createSv;
import static com.hartwig.hmftools.esvee.caller.SvDataCache.buildBreakendMap;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.junit.Test;

public class FiltersTest
{
    @Test
    public void testDeduplication()
    {
        Variant var1 = createSv("01", CHR_1, CHR_2, 100, 200, POS_ORIENT, NEG_ORIENT, "");
        Variant var2 = createSv("01", CHR_1, CHR_2, 100, 200, POS_ORIENT, NEG_ORIENT, "");

        List<Variant> variants = List.of(var1, var2);
        Map<String,List<Breakend>> chrBreakendMap = Maps.newHashMap();
        buildBreakendMap(variants, chrBreakendMap);

        Deduplication.deduplicateVariants(chrBreakendMap);

        assertTrue(var1.isPass());
        assertFalse(var2.isPass());
    }

    /*

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
        Breakend breakend5 = new Breakend(assemblyAlignment, CHR_1, 104, REVERSE, "AAG", homology);
        Breakend breakend6 = new Breakend(assemblyAlignment, CHR_1, 104, FORWARD, "AAG", homology);
        Breakend breakend7 = new Breakend(assemblyAlignment, CHR_1, 104, FORWARD, "AAG", homology);
        Breakend breakend8 = new Breakend(assemblyAlignment, CHR_1, 104, FORWARD, "CCC", homology);

        List<Breakend> breakends = Lists.newArrayList(
                breakend1, breakend2, breakend3, breakend4, breakend5, breakend6, breakend7, breakend8);
        Collections.sort(breakends);

        Deduplication.deduplicateBreakends(breakends);

        assertTrue(breakend1.passing());
        assertTrue(breakend2.passing());
        assertTrue(breakend3.filters().contains(FilterType.DUPLICATE));
        assertFalse(breakend4.filters().contains(FilterType.DUPLICATE));
        assertTrue(breakend5.passing());
        assertTrue(breakend6.passing());
        assertTrue(breakend7.filters().contains(FilterType.DUPLICATE));
        assertTrue(breakend8.passing());
    }
     */

}

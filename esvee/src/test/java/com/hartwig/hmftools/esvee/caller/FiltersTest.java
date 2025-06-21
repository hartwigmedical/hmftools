package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LINKS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.AVG_FRAG_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISC_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.UNIQUE_FRAG_POSITIONS;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.esvee.caller.CallerApplication.isGermline;
import static com.hartwig.hmftools.esvee.caller.CallerTestUtils.TEST_REF_ID;
import static com.hartwig.hmftools.esvee.caller.CallerTestUtils.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.esvee.caller.CallerTestUtils.createSv;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.INV_ADJACENT_MIN_UPS;
import static com.hartwig.hmftools.esvee.caller.SvDataCache.buildBreakendMap;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.common.FilterType;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;
import com.hartwig.hmftools.esvee.prep.types.DiscordantFragType;
import com.hartwig.hmftools.esvee.prep.types.DiscordantStats;

import org.junit.Test;

public class FiltersTest
{
    private static final FilterConstants FILTER_CONSTANTS = FilterConstants.from(false, V37, FilterConstants.DEFAULT_PON_DISTANCE);

    private static final FragmentLengthBounds FRAG_LENGTHS = new FragmentLengthBounds(
            100, 900, 500, 0.1);

    private final VariantFilters mVariantFilters = new VariantFilters(FILTER_CONSTANTS, FRAG_LENGTHS, new DiscordantStats());

    @Test
    public void testDeduplication()
    {
        Variant var1 = createSv("01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "");
        Variant var2 = createSv("02", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "");

        List<Variant> variants = List.of(var1, var2);
        Map<String,List<Breakend>> chrBreakendMap = Maps.newHashMap();
        buildBreakendMap(variants, chrBreakendMap);

        Deduplication.deduplicateVariants(chrBreakendMap);

        assertTrue(var1.isPass());
        assertFalse(var2.isPass());
    }

    @Test
    public void testMinAfFilter()
    {
        Map<String,Object> tumorAttributes = Maps.newHashMap();

        tumorAttributes.put(SPLIT_FRAGS, 1);
        tumorAttributes.put(DISC_FRAGS, 0);
        tumorAttributes.put(TOTAL_FRAGS, 1);
        tumorAttributes.put(REF_DEPTH, 100);

        Map<String,Object> referenceAttributes = Maps.newHashMap();

        referenceAttributes.put(SPLIT_FRAGS, 0); // not above the threshold in the ref sample
        referenceAttributes.put(TOTAL_FRAGS, 0);
        referenceAttributes.put(REF_DEPTH, 100);

        Variant var = createSv(
                "01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "",
                null, referenceAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.MIN_AF));

        tumorAttributes.put(REF_DEPTH, 2000);

        var = createSv(
                "01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "",
                null, referenceAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.MIN_AF));

        // again for a SGL
        tumorAttributes.put(REF_DEPTH, 100);

        var = createSv(
                "01", CHR_1, null, 100, 0, ORIENT_FWD, ORIENT_REV, "",
                null, referenceAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.MIN_AF));

        Map<String,Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(LINE_SITE, true);

        var = createSv(
                "01", CHR_1, null, 100, 0, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, referenceAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.MIN_AF));
    }

    @Test
    public void testMinFragmentLength()
    {
        Map<String,Object> commonAttributes = Maps.newHashMap();

        Variant var = createSv(
                "01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, null, null);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SHORT_FRAG_LENGTH));

        commonAttributes.put(AVG_FRAG_LENGTH, 300);
        commonAttributes.put(TOTAL_FRAGS, 50);

        var = createSv(
                "01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, null, null);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SHORT_FRAG_LENGTH));

        // not tested for LINE or SGLs
        commonAttributes.put(LINE_SITE, true);

        var = createSv(
                "01", CHR_1, CHR_2, 100, 200, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, null, null);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SHORT_FRAG_LENGTH));

        commonAttributes.clear();
        commonAttributes.put(AVG_FRAG_LENGTH, 300);
        commonAttributes.put(TOTAL_FRAGS, 50);

        var = createSv(
                "01", CHR_1, null, 100, 0, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, null, null);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SHORT_FRAG_LENGTH));
    }

    @Test
    public void testInversionShortLowVafHomology()
    {
        Map<String, Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(HOMSEQ, "AGGCT");
        commonAttributes.put(IHOMPOS, new int[] {-10,10});

        Variant var = createSv(
                "01", CHR_1, CHR_1, 100, 200, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, null, null);

        var.contextStart().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 5);
        var.contextEnd().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 5);

        long[] typeCounts = new long[DiscordantFragType.values().length];
        typeCounts[DiscordantFragType.Inv1To5K.ordinal()] = 10;
        DiscordantStats discordantStats = new DiscordantStats(100_000, 10_000, typeCounts);

        VariantFilters invFilters = new VariantFilters(FILTER_CONSTANTS, FRAG_LENGTHS, discordantStats);
        invFilters.applyFilters(var);
        assertFalse(var.filters().contains(FilterType.INV_SHORT_LOW_VAF_HOM));

        typeCounts[DiscordantFragType.Inv1To5K.ordinal()] = 10_000;
        discordantStats = new DiscordantStats(100_000, 10_000, typeCounts);
        invFilters = new VariantFilters(FILTER_CONSTANTS, FRAG_LENGTHS, discordantStats);
        invFilters.applyFilters(var);
        assertTrue(var.filters().contains(FilterType.INV_SHORT_LOW_VAF_HOM));
    }

    @Test
    public void testInversionShortFragmentLowVaf()
    {
        Variant var = createSv(
                "01", CHR_1, CHR_1, 100, 200, ORIENT_FWD, ORIENT_FWD, "",
                Maps.newHashMap(), null, null);

        var.contextStart().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 3);
        var.contextEnd().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 3);

        long[] typeCounts = new long[DiscordantFragType.values().length];
        typeCounts[DiscordantFragType.InvLt1K.ordinal()] = 6_000;
        DiscordantStats discordantStats = new DiscordantStats(100_000, 10_000, typeCounts);

        // set required rate to 2.5%
        VariantFilters invFilters = new VariantFilters(FILTER_CONSTANTS, FRAG_LENGTHS, discordantStats);
        invFilters.applyFilters(var);
        assertTrue(var.filters().contains(FilterType.INV_SHORT_FRAG_LOW_VAF));
    }

    @Test
    public void testShortLowVafDeletion()
    {
        Map<String, Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(HOMSEQ, "AGGCT");
        commonAttributes.put(IHOMPOS, new int[] {-10,10});

        commonAttributes.put(AVG_FRAG_LENGTH, 300);

        Map<String,Object> tumorAttributes = Maps.newHashMap();

        tumorAttributes.put(SPLIT_FRAGS, 5);
        tumorAttributes.put(DISC_FRAGS, 0);
        tumorAttributes.put(TOTAL_FRAGS, 5);
        tumorAttributes.put(REF_DEPTH, 101);

        Variant var = createSv(
                "01", CHR_1, CHR_1, 100, 200, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, null, tumorAttributes);

        var.contextStart().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 5);
        var.contextEnd().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 5);

        mVariantFilters.applyFilters(var);
        assertTrue(var.filters().contains(FilterType.DEL_SHORT_LOW_VAF));
    }

    @Test
    public void testIsolatedInvFilter()
    {
        Map<String, Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(IHOMPOS, new int[] {-10,10});

        Variant inv1 = createSv(
                "01", CHR_1, CHR_1, 100, 150, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, null, Collections.emptyMap());

        // 2 artefacts close to each other, but still marked as isolated
        Variant inv2 = createSv(
                "02", CHR_1, CHR_1, 1100, 1150, ORIENT_REV, ORIENT_REV, "",
                commonAttributes, null, Collections.emptyMap());
        inv2.filters().add(FilterType.PON);

        Variant inv3 = createSv(
                "03", CHR_1, CHR_1, 1140, 1200, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, null, Collections.emptyMap());
        inv3.filters().add(FilterType.INV_SHORT_FRAG_LOW_VAF);

        // chained and higher-UFP INVs are excluded
        commonAttributes = Maps.newHashMap();
        commonAttributes.put(UNIQUE_FRAG_POSITIONS, INV_ADJACENT_MIN_UPS + 1);
        commonAttributes.put(ASM_LINKS, "123");
        Variant inv4 = createSv(
                "04", CHR_1, CHR_1, 2000, 2050, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, null, Collections.emptyMap());

        // lastly adjacent to another breakend
        Variant inv5 = createSv(
                "06", CHR_1, CHR_1, 4000, 4050, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, null, Collections.emptyMap());

        Variant var = createSv(
                "07", CHR_1, CHR_1, 4100, 10000, ORIENT_FWD, ORIENT_REV, "",
                commonAttributes, null, Collections.emptyMap());

        Map<String,List<Breakend>> chrBreakendMap = Maps.newHashMap();
        List<Variant> variants = Lists.newArrayList(inv1, inv2, inv3, inv4, inv5, inv5, var);

        buildBreakendMap(variants, chrBreakendMap);

        mVariantFilters.applyAdjacentFilters(chrBreakendMap);

        assertTrue(inv1.filters().contains(FilterType.INV_SHORT_ISOLATED));
        assertTrue(inv2.filters().contains(FilterType.INV_SHORT_ISOLATED));
        assertTrue(inv3.filters().contains(FilterType.INV_SHORT_ISOLATED));

        assertFalse(inv4.filters().contains(FilterType.INV_SHORT_ISOLATED));
        assertFalse(inv5.filters().contains(FilterType.INV_SHORT_ISOLATED));
        assertFalse(var.filters().contains(FilterType.INV_SHORT_ISOLATED));
    }

    @Test
    public void testMarkGermline()
    {
        Map<String, Object> commonAttributes = Maps.newHashMap();

        Variant var = createSv(
                "01", CHR_1, CHR_1, 100, 200, ORIENT_FWD, ORIENT_FWD, "",
                commonAttributes, null, null);

        var.contextStart().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 50);
        var.contextStart().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(REF_DEPTH, 50);

        var.contextStart().getGenotype(TEST_REF_ID).getExtendedAttributes().put(TOTAL_FRAGS, 2);
        var.contextStart().getGenotype(TEST_REF_ID).getExtendedAttributes().put(REF_DEPTH, 7);

        assertTrue(isGermline(var, TEST_REF_ID));

        var.contextStart().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 500);
        var.contextStart().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(REF_DEPTH, 500);

        var.contextStart().getGenotype(TEST_REF_ID).getExtendedAttributes().put(TOTAL_FRAGS, 2);
        var.contextStart().getGenotype(TEST_REF_ID).getExtendedAttributes().put(REF_DEPTH, 7);

        assertFalse(isGermline(var, TEST_REF_ID));
    }
}

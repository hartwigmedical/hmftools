package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.sv.SvVcfTags.AVG_FRAG_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.DISC_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SPLIT_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.caller.CallerApplication.isGermline;
import static com.hartwig.hmftools.esvee.caller.CallerTestUtils.TEST_REF_ID;
import static com.hartwig.hmftools.esvee.caller.CallerTestUtils.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.esvee.caller.CallerTestUtils.createSv;
import static com.hartwig.hmftools.esvee.caller.SvDataCache.buildBreakendMap;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.esvee.common.FilterType;
import com.hartwig.hmftools.esvee.common.FragmentLengthBounds;

import org.junit.Test;

public class FiltersTest
{
    private static final FilterConstants FILTER_CONSTANTS = FilterConstants.from(false, V37, FilterConstants.DEFAULT_PON_DISTANCE);

    private static final FragmentLengthBounds FRAG_LENGTHS = new FragmentLengthBounds(
            100, 900, 500, 0.1);

    private final VariantFilters mVariantFilters = new VariantFilters(FILTER_CONSTANTS, FRAG_LENGTHS, 0.001);

    @Test
    public void testDeduplication()
    {
        Variant var1 = createSv("01", CHR_1, CHR_2, 100, 200, POS_ORIENT, NEG_ORIENT, "");
        Variant var2 = createSv("02", CHR_1, CHR_2, 100, 200, POS_ORIENT, NEG_ORIENT, "");

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
                "01", CHR_1, CHR_2, 100, 200, POS_ORIENT, NEG_ORIENT, "",
                null, referenceAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.MIN_AF));

        tumorAttributes.put(REF_DEPTH, 2000);

        var = createSv(
                "01", CHR_1, CHR_2, 100, 200, POS_ORIENT, NEG_ORIENT, "",
                null, referenceAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.MIN_AF));

        // again for a SGL
        tumorAttributes.put(REF_DEPTH, 100);

        var = createSv(
                "01", CHR_1, null, 100, 0, POS_ORIENT, NEG_ORIENT, "",
                null, referenceAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.MIN_AF));

        Map<String,Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(LINE_SITE, true);

        var = createSv(
                "01", CHR_1, null, 100, 0, POS_ORIENT, NEG_ORIENT, "",
                commonAttributes, referenceAttributes, tumorAttributes);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.MIN_AF));
    }

    @Test
    public void testMinFragmentLength()
    {
        Map<String,Object> commonAttributes = Maps.newHashMap();

        Variant var = createSv(
                "01", CHR_1, CHR_2, 100, 200, POS_ORIENT, NEG_ORIENT, "",
                commonAttributes, null, null);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SHORT_FRAG_LENGTH));

        commonAttributes.put(AVG_FRAG_LENGTH, 300);
        commonAttributes.put(TOTAL_FRAGS, 50);

        var = createSv(
                "01", CHR_1, CHR_2, 100, 200, POS_ORIENT, NEG_ORIENT, "",
                commonAttributes, null, null);

        mVariantFilters.applyFilters(var);

        assertTrue(var.filters().contains(FilterType.SHORT_FRAG_LENGTH));

        // not tested for LINE or SGLs
        commonAttributes.put(LINE_SITE, true);

        var = createSv(
                "01", CHR_1, CHR_2, 100, 200, POS_ORIENT, NEG_ORIENT, "",
                commonAttributes, null, null);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SHORT_FRAG_LENGTH));

        commonAttributes.clear();
        commonAttributes.put(AVG_FRAG_LENGTH, 300);
        commonAttributes.put(TOTAL_FRAGS, 50);

        var = createSv(
                "01", CHR_1, null, 100, 0, POS_ORIENT, NEG_ORIENT, "",
                commonAttributes, null, null);

        mVariantFilters.applyFilters(var);

        assertFalse(var.filters().contains(FilterType.SHORT_FRAG_LENGTH));
    }

    @Test
    public void testShortLowVafInversion()
    {
        Map<String, Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(HOMSEQ, "AGGCT");
        commonAttributes.put(IHOMPOS, new int[] {-10,10});

        Variant var = createSv(
                "01", CHR_1, CHR_1, 100, 200, POS_ORIENT, POS_ORIENT, "",
                commonAttributes, null, null);

        var.contextStart().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 5);
        var.contextEnd().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 5);

        VariantFilters invFilters = new VariantFilters(FILTER_CONSTANTS, FRAG_LENGTHS, 0.0001);
        invFilters.applyFilters(var);
        assertFalse(var.filters().contains(FilterType.SHORT_LOW_VAF_INV));

        invFilters = new VariantFilters(FILTER_CONSTANTS, FRAG_LENGTHS, 0.1);
        invFilters.applyFilters(var);
        assertTrue(var.filters().contains(FilterType.SHORT_LOW_VAF_INV));
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
                "01", CHR_1, CHR_1, 100, 200, POS_ORIENT, NEG_ORIENT, "",
                commonAttributes, null, tumorAttributes);

        var.contextStart().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 5);
        var.contextEnd().getGenotype(TEST_SAMPLE_ID).getExtendedAttributes().put(TOTAL_FRAGS, 5);

        mVariantFilters.applyFilters(var);
        assertTrue(var.filters().contains(FilterType.SHORT_LOW_VAF_DEL));
    }
    @Test
    public void testMarkGermline()
    {
        Map<String, Object> commonAttributes = Maps.newHashMap();

        Variant var = createSv(
                "01", CHR_1, CHR_1, 100, 200, POS_ORIENT, POS_ORIENT, "",
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

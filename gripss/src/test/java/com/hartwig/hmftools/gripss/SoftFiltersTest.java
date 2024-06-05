package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.region.ExcludedRegions.getPolyGRegion;
import static com.hartwig.hmftools.common.sv.LineElements.POLY_A_HOMOLOGY;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.QUAL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.STRAND_BIAS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_ASRP;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_ASSR;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.GRIDSS_BAQ;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.INDEL_COUNT;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.READ_PAIRS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SPLIT_READS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.SV_FRAG_COUNT;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestApp.TEST_REF_ID;
import static com.hartwig.hmftools.gripss.GripssTestApp.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.gripss.GripssTestUtils.LINE_INSERT_SEQ_T;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSgl;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSv;
import static com.hartwig.hmftools.gripss.GripssTestUtils.defaultFilterConstants;
import static com.hartwig.hmftools.gripss.filters.FilterType.DISCORDANT_PAIR_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.IMPRECISE;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_HOM_LENGTH_SHORT_INV;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_POLY_A_HOM_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_POLY_G_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_NORMAL_COVERAGE;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_QUAL;
import static com.hartwig.hmftools.gripss.filters.FilterType.MIN_TUMOR_AF;
import static com.hartwig.hmftools.gripss.filters.FilterType.MODIFIED_AF;
import static com.hartwig.hmftools.gripss.filters.FilterType.QUAL_PER_AD;
import static com.hartwig.hmftools.gripss.filters.FilterType.SGL_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_DEL_INS_ARTIFACT;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_SR_NORMAL;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_SR_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.SHORT_STRAND_BIAS;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.sv.gridss.GridssVcfTags;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterType;
import com.hartwig.hmftools.gripss.filters.SoftFilters;

import org.junit.Test;

public class SoftFiltersTest
{
    private final SoftFilters mSoftFilters;
    private final FilterCache mFilterCache;

    private final GenotypeIds mGenotypeIds;
    private  final VcfIdGenerator mIdGenerator;

    public SoftFiltersTest()
    {
        mSoftFilters = new SoftFilters(defaultFilterConstants(), false);
        mFilterCache = new FilterCache();
        mIdGenerator = new VcfIdGenerator();
        mGenotypeIds = new GenotypeIds(0, 1, TEST_REF_ID, TEST_SAMPLE_ID);
    }

    private SvData createLongDel(
            final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides, final Map<String,Object> tumorOverrides)
    {
        return createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 10000, POS_ORIENT, NEG_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);
    }

    private SvData createShortDel(
            final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides, final Map<String,Object> tumorOverrides)
    {
        return createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 200, POS_ORIENT, NEG_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);
    }

    private SvData createSingle(
            final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides, final Map<String,Object> tumorOverrides)
    {
        return createSgl(
                mIdGenerator.nextEventId(), CHR_1, 100, POS_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);
    }

    private static void resetOverrides(
            final Map<String,Object> commonOverrides, final Map<String,Object> refOverrides, final Map<String,Object> tumorOverrides)
    {
        commonOverrides.clear();
        refOverrides.clear();
        tumorOverrides.clear();
    }
    
    private void applyFilters(final SvData sv)
    {
        mFilterCache.clear();
        mSoftFilters.applyFilters(sv, mFilterCache);
    }
    
    private boolean hasFilter(final Breakend breakend, final FilterType filter)
    {
        List<FilterType> filters = mFilterCache.getBreakendFilters(breakend);
        return filters != null && filters.contains(filter);
    }

    @Test
    public void testPassingVariant()
    {
        SvData sv = createLongDel(null, null, null);

        // first check that the defaults pass all soft filters
        applyFilters(sv);

        assertNull(mFilterCache.getBreakendFilters(sv.breakendStart()));
        assertNull(mFilterCache.getBreakendFilters(sv.breakendEnd()));
    }

    @Test
    public void testNormalSoftFilters()
    {
        SvData sv = createLongDel(null, null, null);

        Map<String,Object> commonOverrides = Maps.newHashMap();
        Map<String,Object> refOverrides = Maps.newHashMap();
        Map<String,Object> tumorOverrides = Maps.newHashMap();

        // normalCoverageFilter
        refOverrides.put(REF_DEPTH, 1);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MIN_NORMAL_COVERAGE));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // normal relative support
        refOverrides.put(SV_FRAG_COUNT, 80);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MAX_NORMAL_RELATIVE_SUPPORT));



    }


    @Test
    public void testVafSoftFilters()
    {
        SvData sv = createLongDel(null, null, null);

        Map<String,Object> commonOverrides = Maps.newHashMap();
        Map<String,Object> refOverrides = Maps.newHashMap();
        Map<String,Object> tumorOverrides = Maps.newHashMap();

        // min allele frequency
        tumorOverrides.put(SV_FRAG_COUNT, 2);
        tumorOverrides.put(REF_DEPTH, 1000);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MIN_TUMOR_AF));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // min qual
        tumorOverrides.put(QUAL, 1);
        tumorOverrides.put(GRIDSS_BAQ, 1);

        sv = createSingle(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MIN_QUAL));

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MIN_QUAL));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // modified AF
        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        tumorOverrides.put(SV_FRAG_COUNT, 1);
        tumorOverrides.put(INDEL_COUNT, 1);
        tumorOverrides.put(REF_DEPTH, 50);
        tumorOverrides.put(REF_DEPTH_PAIR, 50);

        sv = createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 10000, POS_ORIENT, NEG_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MODIFIED_AF));

        // qual per AD
        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        tumorOverrides.put(SV_FRAG_COUNT, 100);
        tumorOverrides.put(INDEL_COUNT, 100); // thereby lowering qual per tumor fragment

        sv = createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 10000, POS_ORIENT, NEG_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), QUAL_PER_AD));

    }

    @Test
    public void testPositionAndSequenceFilters()
    {
        SvData sv = createLongDel(null, null, null);

        Map<String,Object> commonOverrides = Maps.newHashMap();
        Map<String,Object> refOverrides = Maps.newHashMap();
        Map<String,Object> tumorOverrides = Maps.newHashMap();

        // short split read tumor
        tumorOverrides.put(SPLIT_READS, 0);
        tumorOverrides.put(INDEL_COUNT, 0);
        tumorOverrides.put(GRIDSS_ASSR, 0);

        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), SHORT_SR_SUPPORT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // short split read normal
        refOverrides.put(SPLIT_READS, 1);
        refOverrides.put(INDEL_COUNT, 1);

        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), SHORT_SR_NORMAL));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // discordant support
        refOverrides.put(READ_PAIRS, 0);
        refOverrides.put(GRIDSS_ASRP, 0);
        refOverrides.put(GRIDSS_ASSR, 0);
        tumorOverrides.put(READ_PAIRS, 0);
        tumorOverrides.put(GRIDSS_ASRP, 0);
        tumorOverrides.put(GRIDSS_ASSR, 0);

        sv = createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 148, POS_ORIENT, POS_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), DISCORDANT_PAIR_SUPPORT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // single strand bias
        commonOverrides.put(STRAND_BIAS, 0.99);

        sv = createSingle(commonOverrides, refOverrides, tumorOverrides);

        SvData svLine = createSgl(
                mIdGenerator.nextEventId(), CHR_1, 100, POS_ORIENT, LINE_INSERT_SEQ_T, mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), SGL_STRAND_BIAS));

        applyFilters(svLine);
        assertFalse(hasFilter(svLine.breakendStart(), SGL_STRAND_BIAS));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // single insert sequence min length
        commonOverrides.put(STRAND_BIAS, 0.99);

        sv = createSingle(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), SGL_STRAND_BIAS));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // short DEL artefact
        sv = createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 110, POS_ORIENT, NEG_ORIENT, "AGCTAGCTA",
                mGenotypeIds, commonOverrides, refOverrides, tumorOverrides);

        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), SHORT_DEL_INS_ARTIFACT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // imprecise
        commonOverrides.put(GridssVcfTags.IMPRECISE, "true");

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), IMPRECISE));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // MAX_POLY_G_LENGTH
        ChrBaseRegion excludedRegion = getPolyGRegion(RefGenomeVersion.V37);
        sv = createSv(
                mIdGenerator.nextEventId(), excludedRegion.Chromosome, excludedRegion.Chromosome,
                excludedRegion.start() + 1, excludedRegion.end() - 1, POS_ORIENT, NEG_ORIENT,
                "", mGenotypeIds, commonOverrides, refOverrides, tumorOverrides);

        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MAX_POLY_G_LENGTH));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // MAX_POLY_A_HOM_LENGTH
        commonOverrides.put(HOMSEQ, POLY_A_HOMOLOGY);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MAX_POLY_A_HOM_LENGTH));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // MAX_HOM_LENGTH_SHORT_INV
        commonOverrides.put(HOMSEQ, "AGCGATAA");

        sv = createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 130, POS_ORIENT, POS_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MAX_HOM_LENGTH_SHORT_INV));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // strandBias
        commonOverrides.put(STRAND_BIAS, 0.99);

        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), SHORT_STRAND_BIAS));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // min length
        sv = createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 120, POS_ORIENT, NEG_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MIN_LENGTH));

    }
}

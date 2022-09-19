package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.sv.ExcludedRegions.getPolyGRegion;
import static com.hartwig.hmftools.common.sv.LineElements.POLY_A_HOMOLOGY;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestApp.TEST_REF_ID;
import static com.hartwig.hmftools.gripss.GripssTestApp.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_1;
import static com.hartwig.hmftools.gripss.GripssTestUtils.LINE_INSERT_SEQ_T;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSgl;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSv;
import static com.hartwig.hmftools.gripss.GripssTestUtils.defaultFilterConstants;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_ASRP;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_ASSR;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_BAQ;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_HOMSEQ;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IC;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IMPRECISE;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_QUAL;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REF;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_RP;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_SB;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_SR;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_VF;
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
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
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
        mSoftFilters = new SoftFilters(defaultFilterConstants());
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
    public void testSoftFilters()
    {
        SvData sv = createLongDel(null, null, null);

        // first check that the defaults pass all soft filters
        applyFilters(sv);

        assertNull(mFilterCache.getBreakendFilters(sv.breakendStart()));
        assertNull(mFilterCache.getBreakendFilters(sv.breakendEnd()));

        Map<String,Object> commonOverrides = Maps.newHashMap();
        Map<String,Object> refOverrides = Maps.newHashMap();
        Map<String,Object> tumorOverrides = Maps.newHashMap();

        // normalCoverageFilter
        refOverrides.put(VT_REF, 1);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MIN_NORMAL_COVERAGE));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // normal relative support
        refOverrides.put(VT_VF, 80);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MAX_NORMAL_RELATIVE_SUPPORT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // min allele frequency
        tumorOverrides.put(VT_VF, 2);
        tumorOverrides.put(VT_REF, 1000);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MIN_TUMOR_AF));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // min qual
        tumorOverrides.put(VT_QUAL, 1);
        tumorOverrides.put(VT_BAQ, 1);

        sv = createSingle(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MIN_QUAL));

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MIN_QUAL));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // short split read tumor
        tumorOverrides.put(VT_SR, 0);
        tumorOverrides.put(VT_IC, 0);
        tumorOverrides.put(VT_ASSR, 0);

        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), SHORT_SR_SUPPORT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // short split read normal
        refOverrides.put(VT_SR, 1);
        refOverrides.put(VT_IC, 1);

        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), SHORT_SR_NORMAL));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // discordant support
        refOverrides.put(VT_RP, 0);
        refOverrides.put(VT_ASRP, 0);
        refOverrides.put(VT_ASSR, 0);
        tumorOverrides.put(VT_RP, 0);
        tumorOverrides.put(VT_ASRP, 0);
        tumorOverrides.put(VT_ASSR, 0);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), DISCORDANT_PAIR_SUPPORT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // single strand bias
        commonOverrides.put(VT_SB, 0.99);

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
        commonOverrides.put(VT_SB, 0.99);

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
        commonOverrides.put(VT_IMPRECISE, "true");

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
        commonOverrides.put(VT_HOMSEQ, POLY_A_HOMOLOGY);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MAX_POLY_A_HOM_LENGTH));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // MAX_HOM_LENGTH_SHORT_INV
        commonOverrides.put(VT_HOMSEQ, "AGCGATAA");

        sv = createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 130, POS_ORIENT, POS_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        applyFilters(sv);
        assertTrue(hasFilter(sv.breakendStart(), MAX_HOM_LENGTH_SHORT_INV));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // strandBias
        commonOverrides.put(VT_SB, 0.99);

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

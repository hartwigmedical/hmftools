package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_REF_ID;
import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_1;
import static com.hartwig.hmftools.gripss.GripssTestUtils.LINE_INSERT_SEQ_T;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSgl;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSv;
import static com.hartwig.hmftools.gripss.GripssTestUtils.defaultFilterConstants;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_ASRP;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_BAQ;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_BASRP;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_BASSR;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_HOMSEQ;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IC;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IHOMPOS;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IMPRECISE;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_QUAL;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REF;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REFPAIR;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_SB;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_SR;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_VF;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.LINC_00486_V37;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.POLY_A_HOMOLOGY;
import static com.hartwig.hmftools.gripss.filters.FilterType.DISCORDANT_PAIR_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterType.IMPRECISE;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_HOM_LENGTH_SHORT_INV;
import static com.hartwig.hmftools.gripss.filters.FilterType.MAX_INEXACT_HOM_LENGTH_SHORT_DEL;
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
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.SoftFilters;

import org.junit.Test;

public class SoftFiltersTest
{
    private final SoftFilters mSoftFilters;
    private final GenotypeIds mGenotypeIds;
    private  final VcfIdGenerator mIdGenerator;

    public SoftFiltersTest()
    {
        mSoftFilters = new SoftFilters(defaultFilterConstants());
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

    @Test
    public void testSoftFilters()
    {
        SvData sv = createLongDel(null, null, null);

        // first check that the defaults pass all soft filters
        mSoftFilters.applyFilters(sv);

        assertTrue(sv.breakendStart().getFilters().isEmpty());
        assertTrue(sv.getFilters().isEmpty());

        Map<String,Object> commonOverrides = Maps.newHashMap();
        Map<String,Object> refOverrides = Maps.newHashMap();
        Map<String,Object> tumorOverrides = Maps.newHashMap();

        // normalCoverageFilter
        refOverrides.put(VT_REF, 1);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(MIN_NORMAL_COVERAGE));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // normal relative support
        refOverrides.put(VT_VF, 80);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(MAX_NORMAL_RELATIVE_SUPPORT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // min allele frequency
        tumorOverrides.put(VT_VF, 2);
        refOverrides.put(VT_REF, 1000);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(MIN_TUMOR_AF));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // min qual
        tumorOverrides.put(VT_QUAL, 1);
        tumorOverrides.put(VT_BAQ, 1);

        sv = createSingle(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(MIN_QUAL));

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(MIN_QUAL));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // short split read tumor
        tumorOverrides.put(VT_SR, 0);
        tumorOverrides.put(VT_IC, 0);

        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(SHORT_SR_SUPPORT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // short split read normal
        refOverrides.put(VT_SR, 1);
        refOverrides.put(VT_IC, 1);

        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(SHORT_SR_NORMAL));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // discordant support
        refOverrides.put(VT_REFPAIR, 0);
        refOverrides.put(VT_ASRP, 0);
        tumorOverrides.put(VT_REFPAIR, 0);
        tumorOverrides.put(VT_ASRP, 0);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(DISCORDANT_PAIR_SUPPORT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // single strand bias
        commonOverrides.put(VT_SB, 0.99);

        sv = createSingle(commonOverrides, refOverrides, tumorOverrides);

        SvData svLine = createSgl(
                mIdGenerator.nextEventId(), CHR_1, 100, POS_ORIENT, LINE_INSERT_SEQ_T, mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(SGL_STRAND_BIAS));

        mSoftFilters.applyFilters(svLine);
        assertFalse(svLine.breakendStart().getFilters().contains(SGL_STRAND_BIAS));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // single insert sequence min length
        commonOverrides.put(VT_SB, 0.99);

        sv = createSingle(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(SGL_STRAND_BIAS));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // short DEL artefact
        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(SHORT_DEL_INS_ARTIFACT));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // imprecise
        commonOverrides.put(VT_IMPRECISE, "true");

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.getFilters().contains(IMPRECISE));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // MAX_POLY_G_LENGTH
        sv = createSv(
                mIdGenerator.nextEventId(), LINC_00486_V37.Chromosome, LINC_00486_V37.Chromosome,
                LINC_00486_V37.start() + 1, LINC_00486_V37.end() - 1, POS_ORIENT, NEG_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        mSoftFilters.applyFilters(sv);
        assertTrue(sv.getFilters().contains(MAX_POLY_G_LENGTH));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // MAX_POLY_A_HOM_LENGTH
        commonOverrides.put(VT_HOMSEQ, POLY_A_HOMOLOGY);

        sv = createLongDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.getFilters().contains(MAX_POLY_A_HOM_LENGTH));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // MAX_HOM_LENGTH_SHORT_INV
        commonOverrides.put(VT_HOMSEQ, "AGCGATAA");

        sv = createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 130, POS_ORIENT, POS_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        mSoftFilters.applyFilters(sv);
        assertTrue(sv.getFilters().contains(MAX_HOM_LENGTH_SHORT_INV));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // inexactHomologyLengthShortDel
        commonOverrides.put(VT_IHOMPOS, Lists.newArrayList(-5, 8));

        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(MAX_INEXACT_HOM_LENGTH_SHORT_DEL));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // strandBias
        commonOverrides.put(VT_SB, 0.99);

        sv = createShortDel(commonOverrides, refOverrides, tumorOverrides);
        mSoftFilters.applyFilters(sv);
        assertTrue(sv.breakendStart().getFilters().contains(SHORT_STRAND_BIAS));

        resetOverrides(commonOverrides, refOverrides, tumorOverrides);

        // min length
        sv = createSv(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 120, POS_ORIENT, NEG_ORIENT, "", mGenotypeIds,
                commonOverrides, refOverrides, tumorOverrides);

        mSoftFilters.applyFilters(sv);
        assertTrue(sv.getFilters().contains(MIN_LENGTH));
    }




}

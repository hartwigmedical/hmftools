package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_REF_ID;
import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_1;
import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_2;
import static com.hartwig.hmftools.gripss.GripssTestUtils.DEFAULT_QUAL;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSglBreakend;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSvBreakends;
import static com.hartwig.hmftools.gripss.GripssTestUtils.defaultFilterConstants;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_QUAL;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertNull;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.HotspotCache;
import com.hartwig.hmftools.gripss.filters.KnownHotspot;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantBuilderTest
{
    private final HotspotCache mHotspotCache;
    private final VariantBuilder mBuilder;
    private  final VcfIdGenerator mIdGenerator;
    private final GenotypeIds mGenotypeIds;

    public VariantBuilderTest()
    {
        mHotspotCache = new HotspotCache(null);
        mBuilder = new VariantBuilder(defaultFilterConstants(), mHotspotCache);
        mIdGenerator = new VcfIdGenerator();
        mGenotypeIds = new GenotypeIds(0, 1, TEST_REF_ID, TEST_SAMPLE_ID);
    }

    @Test
    public void testVariantTypeCreation()
    {
        // each SV in turn
        VariantContext[] contexts = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 200, POS_ORIENT, NEG_ORIENT, "A", "");

        SvData del = mBuilder.checkCreateVariant(contexts[SE_START], mGenotypeIds);
        assertNull(del);
        del = mBuilder.checkCreateVariant(contexts[SE_END], mGenotypeIds);
        assertNotNull(del);

        assertEquals(CHR_1, del.chromosomeStart());
        assertEquals(CHR_1, del.chromosomeEnd());
        assertEquals(100, del.posStart());
        assertEquals(200, del.posEnd());
        assertEquals(POS_ORIENT, del.orientStart());
        assertEquals(NEG_ORIENT, del.orientEnd());
        assertEquals(DEL, del.type());

        contexts = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 200, NEG_ORIENT, POS_ORIENT, "A", "");

        SvData dup = mBuilder.checkCreateVariant(contexts[SE_START], mGenotypeIds);
        assertNull(dup);
        dup = mBuilder.checkCreateVariant(contexts[SE_END], mGenotypeIds);
        assertNotNull(dup);

        assertEquals(CHR_1, dup.chromosomeStart());
        assertEquals(CHR_1, dup.chromosomeEnd());
        assertEquals(100, dup.posStart());
        assertEquals(200, dup.posEnd());
        assertEquals(NEG_ORIENT, dup.orientStart());
        assertEquals(POS_ORIENT, dup.orientEnd());
        assertEquals(DUP, dup.type());

        contexts = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 200, NEG_ORIENT, NEG_ORIENT, "A", "");

        SvData inv = mBuilder.checkCreateVariant(contexts[SE_START], mGenotypeIds);
        assertNull(inv);
        inv = mBuilder.checkCreateVariant(contexts[SE_END], mGenotypeIds);
        assertNotNull(inv);

        assertEquals(CHR_1, inv.chromosomeStart());
        assertEquals(CHR_1, inv.chromosomeEnd());
        assertEquals(100, inv.posStart());
        assertEquals(200, inv.posEnd());
        assertEquals(NEG_ORIENT, inv.orientStart());
        assertEquals(NEG_ORIENT, inv.orientEnd());
        assertEquals(INV, inv.type());

        contexts = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_2, 100, 200, NEG_ORIENT, NEG_ORIENT, "A", "");

        SvData bnd = mBuilder.checkCreateVariant(contexts[SE_START], mGenotypeIds);
        assertNull(bnd);
        bnd = mBuilder.checkCreateVariant(contexts[SE_END], mGenotypeIds);
        assertNotNull(bnd);

        assertEquals(CHR_1, bnd.chromosomeStart());
        assertEquals(CHR_2, bnd.chromosomeEnd());
        assertEquals(100, bnd.posStart());
        assertEquals(200, bnd.posEnd());
        assertEquals(NEG_ORIENT, bnd.orientStart());
        assertEquals(NEG_ORIENT, bnd.orientEnd());
        assertEquals(BND, bnd.type());

        VariantContext sglContext = createSglBreakend(mIdGenerator.nextEventId(), CHR_1, 100, POS_ORIENT, "A", "");

        SvData sgl = mBuilder.checkCreateVariant(sglContext, mGenotypeIds);
        assertNotNull(sgl);

        assertEquals(CHR_1, sgl.chromosomeStart());
        assertEquals(100, sgl.posStart());
        assertEquals(POS_ORIENT, sgl.orientStart());
        assertEquals(SGL, sgl.type());

        assertEquals(0, mBuilder.hardFilteredCount());
        assertEquals(0, mBuilder.incompleteSVs());
    }

    @Test
    public void testHardFilters()
    {
        // variant with first leg hard-filtered
        VariantContext[] contexts = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 200, POS_ORIENT, NEG_ORIENT, "A", "");

        contexts[SE_START].getGenotype(1).getExtendedAttributes().put(VT_QUAL, 1);

        SvData sv = mBuilder.checkCreateVariant(contexts[SE_START], mGenotypeIds);
        assertNull(sv);
        sv = mBuilder.checkCreateVariant(contexts[SE_END], mGenotypeIds);
        assertNull(sv);

        assertEquals(1, mBuilder.hardFilteredCount());
        assertEquals(0, mBuilder.incompleteSVs());

        mBuilder.clearState();

        // second leg hard-filtered
        contexts[SE_END].getGenotype(1).getExtendedAttributes().put(VT_QUAL, 1);
        contexts[SE_START].getGenotype(1).getExtendedAttributes().put(VT_QUAL, DEFAULT_QUAL);

        sv = mBuilder.checkCreateVariant(contexts[SE_START], mGenotypeIds);
        assertNull(sv);
        sv = mBuilder.checkCreateVariant(contexts[SE_END], mGenotypeIds);
        assertNull(sv);

        assertEquals(1, mBuilder.hardFilteredCount());
        assertEquals(0, mBuilder.incompleteSVs());

        mBuilder.clearState();

        // SGL hard-filtered
        VariantContext sglContext = createSglBreakend(mIdGenerator.nextEventId(), CHR_1, 100, POS_ORIENT, "A", "");
        sglContext.getGenotype(1).getExtendedAttributes().put(VT_QUAL, 1);

        SvData sgl = mBuilder.checkCreateVariant(sglContext, mGenotypeIds);
        assertNull(sgl);

        assertEquals(1, mBuilder.hardFilteredCount());
        assertEquals(0, mBuilder.incompleteSVs());
    }

    @Test
    public void testHotSpotVariantCreation()
    {
        mHotspotCache.addHotspot(new KnownHotspot(
                new ChrBaseRegion(CHR_1, 100, 200), POS_ORIENT,
                new ChrBaseRegion(CHR_1, 1999, 2000), NEG_ORIENT, "geneInfo"));

        // first a low-qual pair but matching the hotspot
        VariantContext[] contexts = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 2000, POS_ORIENT, NEG_ORIENT, "A", "");

        contexts[SE_START].getGenotype(1).getExtendedAttributes().put(VT_QUAL, 1);

        SvData sv = mBuilder.checkCreateVariant(contexts[SE_START], mGenotypeIds);
        assertNull(sv);
        sv = mBuilder.checkCreateVariant(contexts[SE_END], mGenotypeIds);
        assertNotNull(sv);

        assertEquals(0, mBuilder.hardFilteredCount());
        assertEquals(0, mBuilder.incompleteSVs());

        mBuilder.clearState();

        // now with only the first leg matching the hotspot
        contexts = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 3000, POS_ORIENT, NEG_ORIENT, "A", "");

        contexts[SE_START].getGenotype(1).getExtendedAttributes().put(VT_QUAL, 1);

        sv = mBuilder.checkCreateVariant(contexts[SE_START], mGenotypeIds);
        assertNull(sv);
        sv = mBuilder.checkCreateVariant(contexts[SE_END], mGenotypeIds);
        assertNull(sv);

        assertEquals(1, mBuilder.hardFilteredCount());
        assertEquals(0, mBuilder.incompleteSVs());

        mBuilder.clearState();

        // second leg hard-filtered but matching a hotspot
        // now with only the first leg matching the hotspot
        contexts = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 2000, POS_ORIENT, NEG_ORIENT, "A", "");

        contexts[SE_END].getGenotype(1).getExtendedAttributes().put(VT_QUAL, 1);

        sv = mBuilder.checkCreateVariant(contexts[SE_START], mGenotypeIds);
        assertNull(sv);
        sv = mBuilder.checkCreateVariant(contexts[SE_END], mGenotypeIds);
        assertNotNull(sv);

        assertEquals(0, mBuilder.hardFilteredCount());
        assertEquals(0, mBuilder.incompleteSVs());
    }

}

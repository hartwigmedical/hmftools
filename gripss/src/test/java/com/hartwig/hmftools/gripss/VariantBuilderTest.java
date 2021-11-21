package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_REF_ID;
import static com.hartwig.hmftools.gripss.GripssTestApplication.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_1;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSvBreakends;
import static com.hartwig.hmftools.gripss.GripssTestUtils.defaultFilterConstants;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertNull;
import static junit.framework.TestCase.assertTrue;

import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class VariantBuilderTest
{
    // private GripssTestApplication mGripss;
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
        // each leg of an SV
        VariantContext[] sv1 = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 200, POS_ORIENT, NEG_ORIENT, "A", "");

        SvData sv = mBuilder.checkCreateVariant(sv1[SE_START], mGenotypeIds);
        assertNull(sv);
        sv = mBuilder.checkCreateVariant(sv1[SE_END], mGenotypeIds);
        assertNotNull(sv);



    }

    @Test
    public void testHotSpotVariantCreation()
    {
        // each leg of an SV
        VariantContext[] sv1 = createSvBreakends(
                mIdGenerator.nextEventId(), CHR_1, CHR_1, 100, 200, POS_ORIENT, NEG_ORIENT, "A", "");

        SvData sv = mBuilder.checkCreateVariant(sv1[SE_START], mGenotypeIds);
        assertFalse(sv != null);
        sv = mBuilder.checkCreateVariant(sv1[SE_END], mGenotypeIds);
        assertTrue(sv != null);
    }

    @Test
    public void testHardFilters()
    {
        //VariantContext[] sv1 = createSvBreakends(
        //        mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 100, 200, POS_ORIENT, NEG_ORIENT, "A", "");

        // mGripss.processVariant(sv1[SE_START]);
        // mGripss.processVariant(sv1[SE_END]);

        // VariantContext var1 = createBreakend(mGripss.IdGen.next(true), CHR_1, 100, "A", "G");

        // mGripss.processVariant(var1);
    }

}

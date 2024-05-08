package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.CIRPOS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.IMPRECISE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSgl;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSv;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;

import org.junit.Test;

public class RealignmentTest
{
    private final GripssTestApp mGripss;
    private final BreakendRealigner mRealigner;

    public RealignmentTest()
    {
        mGripss = new GripssTestApp();
        mRealigner = new BreakendRealigner(mGripss.RefGenome);
        mGripss.RefGenome.RefGenomeMap.put(CHR_1, generateRandomBases(100));
        mGripss.RefGenome.RefGenomeMap.put(CHR_2, generateRandomBases(100));
    }

    @Test
    public void testVariantRealignment()
    {
        Map<String, Object> attributesStart = Maps.newHashMap();
        attributesStart.put(CIPOS, new int[] {-10, 50});
        attributesStart.put(CIRPOS, new int[] {-20, 40});
        attributesStart.put(IHOMPOS, new int[] {-15, 25});

        Map<String, Object> attributesEnd = Maps.newHashMap();
        attributesEnd.put(CIRPOS, new int[] {-10, 50});
        attributesEnd.put(CIPOS, new int[] {-20, 40});
        attributesEnd.put(IHOMPOS, new int[] {-15, 5});

        SvData var = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_2, 20, 80, POS_ORIENT, NEG_ORIENT, "",
                mGripss.GenotypeIds, attributesStart, attributesEnd);

        Breakend newBreakendStart = mRealigner.realign(var.breakendStart(), var.isSgl(), var.imprecise());
        assertTrue(newBreakendStart.realigned());
        assertEquals(-30, newBreakendStart.ConfidenceInterval.Start);
        assertEquals(30, newBreakendStart.ConfidenceInterval.End);
        assertEquals(40, newBreakendStart.Position);
        assertEquals(-35, newBreakendStart.InexactHomology.Start);
        assertEquals(5, newBreakendStart.InexactHomology.End);

        Breakend newBreakendEnd = mRealigner.realignRemote(var.breakendEnd(), newBreakendStart);
        assertTrue(newBreakendEnd.realigned());
        assertEquals(-30, newBreakendEnd.RemoteConfidenceInterval.Start);
        assertEquals(30, newBreakendEnd.RemoteConfidenceInterval.End);
        assertEquals(80, newBreakendEnd.Position);

        // no realignment for inserts
        var = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_2, 50, 50, POS_ORIENT, NEG_ORIENT, "ATGC",
                mGripss.GenotypeIds, attributesStart, attributesEnd);

        assertFalse(mRealigner.realign(var.breakendStart(), var.isSgl(), var.imprecise()).realigned());

        // single realignment
        attributesStart = Maps.newHashMap();
        attributesStart.put(CIPOS, new int[] {-10, 50});

        SvData sgl = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 50, POS_ORIENT, "",
                mGripss.GenotypeIds, attributesStart, null, null);

        assertFalse(mRealigner.realign(sgl.breakendStart(), sgl.isSgl(), sgl.imprecise()).realigned());

        attributesStart.put(IMPRECISE, "true");

        sgl = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 50, POS_ORIENT, "",
                mGripss.GenotypeIds, attributesStart, null, null);

        Breakend newBreakend = mRealigner.realign(sgl.breakendStart(), sgl.isSgl(), sgl.imprecise());
        assertTrue(newBreakend.realigned());
    }

}

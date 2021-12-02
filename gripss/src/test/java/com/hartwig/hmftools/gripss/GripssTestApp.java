package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.Map;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;

public class GripssTestApp
{
    public final MockRefGenome RefGenome;
    public final VcfIdGenerator IdGen;
    public final GenotypeIds GenotypeIds;

    public static final String TEST_SAMPLE_ID = "TUMOR_ID";
    public static final String TEST_REF_ID = "REF_ID";

    public GripssTestApp()
    {
        IdGen = new VcfIdGenerator();
        RefGenome = new MockRefGenome();
        GenotypeIds = new GenotypeIds(0, 1, TEST_REF_ID, TEST_SAMPLE_ID);
    }

    // convenience builders for each type
    public SvData createDel(
            final String chromosome, int posStart, int posEnd,
            final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        return GripssTestUtils.createSv(
                IdGen.nextEventId(), chromosome, chromosome, posStart, posEnd, POS_ORIENT, NEG_ORIENT, "", GenotypeIds,
                attributesStart, attributesEnd);
    }

    public SvData createDup(
            final String chromosome, int posStart, int posEnd,
            final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        return GripssTestUtils.createSv(
                IdGen.nextEventId(), chromosome, chromosome, posStart, posEnd, NEG_ORIENT, POS_ORIENT, "", GenotypeIds,
                attributesStart, attributesEnd);
    }

    public SvData createInv(
            final String chromosome, int posStart, int posEnd, byte orientation,
            final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        return GripssTestUtils.createSv(
                IdGen.nextEventId(), chromosome, chromosome, posStart, posEnd, orientation, orientation, "", GenotypeIds,
                attributesStart, attributesEnd);
    }

    public SvData createBnd(
            final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        return GripssTestUtils.createSv(
                IdGen.nextEventId(), chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, "", GenotypeIds,
                attributesStart, attributesEnd);
    }
}

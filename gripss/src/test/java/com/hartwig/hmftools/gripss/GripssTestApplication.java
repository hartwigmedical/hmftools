package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_1;
import static com.hartwig.hmftools.gripss.GripssTestUtils.defaultFilterConstants;

import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.filters.FilterConstants;

import htsjdk.variant.variantcontext.VariantContext;

public class GripssTestApplication
{
    public final MockRefGenome RefGenome;
    public final VcfIdGenerator IdGen;
    public final GripssApplication mGripss;
    public final GenotypeIds mGenotypeIds;

    public static final String TEST_SAMPLE_ID = "TUMOR_ID";
    public static final String TEST_REF_ID = "REF_ID";

    public GripssTestApplication()
    {
        IdGen = new VcfIdGenerator();
        RefGenome = new MockRefGenome();

        mGenotypeIds = new GenotypeIds(0, 1, TEST_REF_ID, TEST_SAMPLE_ID);

        GripssConfig config = new GripssConfig(TEST_SAMPLE_ID, TEST_REF_ID, V37, "vcf_file");

        FilterConstants filterConstants = defaultFilterConstants();

        mGripss = new GripssApplication(config, filterConstants, RefGenome, null);
    }

    public void processVariant(final VariantContext variant)
    {
        mGripss.processVariant(variant, mGenotypeIds);
    }

    // convenience builders for each type
    public SvData createDel(
            final String chromosome, int posStart, int posEnd,
            final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        return GripssTestUtils.createSv(
                IdGen.nextEventId(), chromosome, chromosome, posStart, posEnd, POS_ORIENT, NEG_ORIENT, "", mGenotypeIds,
                attributesStart, attributesEnd);
    }

    public SvData createDup(
            final String chromosome, int posStart, int posEnd,
            final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        return GripssTestUtils.createSv(
                IdGen.nextEventId(), chromosome, chromosome, posStart, posEnd, NEG_ORIENT, POS_ORIENT, "", mGenotypeIds,
                attributesStart, attributesEnd);
    }

    public SvData createInv(
            final String chromosome, int posStart, int posEnd, byte orientation,
            final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        return GripssTestUtils.createSv(
                IdGen.nextEventId(), chromosome, chromosome, posStart, posEnd, orientation, orientation, "", mGenotypeIds,
                attributesStart, attributesEnd);
    }

    public SvData createBnd(
            final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final Map<String,Object> attributesStart, final Map<String,Object> attributesEnd)
    {
        return GripssTestUtils.createSv(
                IdGen.nextEventId(), chrStart, chrEnd, posStart, posEnd, orientStart, orientEnd, "", mGenotypeIds,
                attributesStart, attributesEnd);
    }
}

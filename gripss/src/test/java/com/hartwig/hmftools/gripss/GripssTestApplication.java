package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.gripss.GripssTestUtils.defaultFilterConstants;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
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
}

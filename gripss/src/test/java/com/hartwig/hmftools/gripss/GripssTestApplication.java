package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_HARD_MAX_NORMAL_ABSOLUTE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_HARD_MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MAX_HOM_LENGTH_SHORT_INV;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MAX_INEXACT_HOM_LENGTH_SHORT_DEL;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MAX_SHORT_STRAND_BIAS;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_LENGTH;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_NORMAL_COVERAGE;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_QUAL_BREAK_END;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_QUAL_BREAK_POINT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_MIN_TUMOR_AF;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_PON_DISTANCE;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.DEFAULT_SOFT_MAX_NORMAL_RELATIVE_SUPPORT;
import static com.hartwig.hmftools.gripss.filters.FilterConstants.LINC_00486_V37;

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

        FilterConstants filterConstants = new FilterConstants(
                FilterConstants.DEFAULT_HARD_MIN_TUMOR_QUAL,
                DEFAULT_HARD_MAX_NORMAL_ABSOLUTE_SUPPORT,
                DEFAULT_HARD_MAX_NORMAL_RELATIVE_SUPPORT,
                DEFAULT_SOFT_MAX_NORMAL_RELATIVE_SUPPORT,
                DEFAULT_MIN_NORMAL_COVERAGE,
                DEFAULT_MIN_TUMOR_AF,
                DEFAULT_MAX_SHORT_STRAND_BIAS,
                DEFAULT_MIN_QUAL_BREAK_END,
                DEFAULT_MIN_QUAL_BREAK_POINT,
                DEFAULT_MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION,
                DEFAULT_MAX_HOM_LENGTH_SHORT_INV,
                DEFAULT_MAX_INEXACT_HOM_LENGTH_SHORT_DEL,
                DEFAULT_MIN_LENGTH,
                DEFAULT_PON_DISTANCE,
                LINC_00486_V37);

        mGripss = new GripssApplication(config, filterConstants, RefGenome, null);
    }

    public void processVariant(final VariantContext variant)
    {
        mGripss.processVariant(variant, mGenotypeIds);
    }
}

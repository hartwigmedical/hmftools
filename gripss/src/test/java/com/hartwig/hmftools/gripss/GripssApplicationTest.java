package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.gripss.GripssTestApp.TEST_REF_ID;
import static com.hartwig.hmftools.gripss.GripssTestApp.TEST_SAMPLE_ID;
import static com.hartwig.hmftools.gripss.GripssTestUtils.defaultFilterConstants;

import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.gripss.filters.FilterConstants;

public class GripssApplicationTest
{
    public final GripssApplication mGripss;
    public final MockRefGenome RefGenome;

    public GripssApplicationTest()
    {
        GripssConfig config = new GripssConfig(TEST_SAMPLE_ID, TEST_REF_ID, V37, "vcf_file");
        FilterConstants filterConstants = defaultFilterConstants();
        RefGenome = new MockRefGenome();

        mGripss = new GripssApplication(config, filterConstants, RefGenome, null);

    }
}

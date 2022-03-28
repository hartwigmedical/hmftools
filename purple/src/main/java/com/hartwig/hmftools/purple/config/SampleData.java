package com.hartwig.hmftools.purple.config;

import com.hartwig.hmftools.purple.SomaticVariantCache;
import com.hartwig.hmftools.purple.StructuralVariantCache;

public class SampleData
{
    public final String ReferenceId;
    public final String SampleId;

    public final AmberData Amber;
    public final CobaltData Cobalt;
    public final StructuralVariantCache SvCache;
    public final SomaticVariantCache SomaticCache;

    public SampleData(
            final String referenceId, final String sampleId, final AmberData amber, final CobaltData cobalt,
            final StructuralVariantCache svCache, final SomaticVariantCache somaticCache)
    {
        ReferenceId = referenceId;
        SampleId = sampleId;

        Amber = amber;
        Cobalt = cobalt;
        SvCache = svCache;
        SomaticCache = somaticCache;
    }
}

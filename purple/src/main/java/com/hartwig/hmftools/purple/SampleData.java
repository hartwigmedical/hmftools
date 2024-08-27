package com.hartwig.hmftools.purple;

import com.hartwig.hmftools.purple.somatic.SomaticVariantCache;
import com.hartwig.hmftools.purple.sv.SomaticSvCache;

public class SampleData
{
    public final String ReferenceId;
    public final String SampleId;

    public final AmberData Amber;
    public final CobaltData Cobalt;
    public final SomaticSvCache SvCache;
    public final SomaticVariantCache SomaticCache;

    public SampleData(
            final String referenceId, final String sampleId, final AmberData amber, final CobaltData cobalt,
            final SomaticSvCache svCache, final SomaticVariantCache somaticCache)
    {
        ReferenceId = referenceId;
        SampleId = sampleId;

        Amber = amber;
        Cobalt = cobalt;
        SvCache = svCache;
        SomaticCache = somaticCache;
    }
}

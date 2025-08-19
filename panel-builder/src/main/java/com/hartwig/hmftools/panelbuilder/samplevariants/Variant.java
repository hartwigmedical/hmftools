package com.hartwig.hmftools.panelbuilder.samplevariants;

import com.hartwig.hmftools.panelbuilder.ProbeTarget;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;

public interface Variant
{
    ProbeTarget generateProbeTarget();

    boolean isDriver();

    TargetMetadata.Type targetType();
}

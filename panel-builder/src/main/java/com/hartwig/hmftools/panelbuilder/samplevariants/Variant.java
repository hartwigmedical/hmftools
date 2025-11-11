package com.hartwig.hmftools.panelbuilder.samplevariants;

import com.hartwig.hmftools.panelbuilder.SequenceDefinition;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;

public interface Variant
{
    SequenceDefinition generateProbe();

    boolean isDriver();

    TargetMetadata.Type targetType();

    String toString();
}

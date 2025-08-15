package com.hartwig.hmftools.panelbuilder.samplevariants;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;

public interface Variant
{
    VariantProbeData generateProbe(final RefGenomeInterface refGenome);

    boolean isDriver();

    TargetMetadata.Type targetType();
}

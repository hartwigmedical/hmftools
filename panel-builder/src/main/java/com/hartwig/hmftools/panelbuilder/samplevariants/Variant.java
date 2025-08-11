package com.hartwig.hmftools.panelbuilder.samplevariants;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;

public interface Variant
{
    VariantProbeData generateProbe(final RefGenomeInterface refGenome);

    boolean isDriver();

    // TODO: is this actually needed?
    List<ProximateLocations.Location> checkedLocations();

    TargetMetadata.Type targetType();
}

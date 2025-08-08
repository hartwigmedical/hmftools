package com.hartwig.hmftools.panelbuilder.samplevariants;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public abstract class Variant
{
    public abstract VariantProbeData generateProbe(final RefGenomeInterface refGenome);

    public abstract boolean isDriver();

    public abstract List<ProximateLocations.Location> checkedLocations();
}

package com.hartwig.hmftools.common.gene;

import java.util.Collection;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

import org.jetbrains.annotations.NotNull;

public class GeneModel {

    @NotNull
    private final Collection<HmfGenomeRegion> regions;
    @NotNull
    private final Set<String> panel;

    public GeneModel(@NotNull Collection<HmfGenomeRegion> regions) {
        this.regions = regions;
        this.panel = regions.stream().map(TranscriptRegion::gene).collect(Collectors.toSet());
    }

    @NotNull
    public Collection<HmfGenomeRegion> regions() {
        return regions;
    }

    public long numberOfBases() {
        return regions.stream().mapToLong(GenomeRegion::bases).sum();
    }

    public int numberOfRegions() {
        return regions.size();
    }

    @NotNull
    public Set<String> panel() {
        return panel;
    }
}

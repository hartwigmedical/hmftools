package com.hartwig.hmftools.patientreporter.algo;

import java.util.Collection;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.region.TranscriptRegion;

import org.jetbrains.annotations.NotNull;

public class GeneModel {

    @NotNull
    private final Collection<HmfTranscriptRegion> regions;
    @NotNull
    private final Set<String> panel;

    public GeneModel(@NotNull Collection<HmfTranscriptRegion> regions) {
        this.regions = regions;
        this.panel = regions.stream().map(TranscriptRegion::gene).collect(Collectors.toSet());
    }

    @NotNull
    public Collection<HmfTranscriptRegion> regions() {
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

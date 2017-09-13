package com.hartwig.hmftools.common.gene;

import java.util.Collection;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.collect.SortedSetMultimap;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;

import org.jetbrains.annotations.NotNull;

public class GeneModel {

    private final Map<String, HmfGenomeRegion> transcriptMap;
    private final Collection<HmfGenomeRegion> regions;
    private final Slicer slicer;

    public GeneModel(@NotNull final SortedSetMultimap<String, HmfGenomeRegion> regions) {
        this.regions = regions.values();
        slicer = SlicerFactory.fromRegions(regions);
        transcriptMap = extractTranscriptMap(regions.values());
    }

    @NotNull
    public Collection<HmfGenomeRegion> hmfRegions() {
        return regions;
    }

    @NotNull
    public Slicer slicer() {
        return slicer;
    }

    @NotNull
    public Map<String, HmfGenomeRegion> transcriptMap() {
        return transcriptMap;
    }

    @NotNull
    private static Map<String, HmfGenomeRegion> extractTranscriptMap(final @NotNull Collection<HmfGenomeRegion> regions) {
        final Map<String, HmfGenomeRegion> transcriptMap = Maps.newHashMap();
        for (final HmfGenomeRegion region : regions) {
            transcriptMap.put(region.transcriptID(), region);
        }
        return transcriptMap;
    }
}

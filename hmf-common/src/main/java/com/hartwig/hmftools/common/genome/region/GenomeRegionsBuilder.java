package com.hartwig.hmftools.common.genome.region;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class GenomeRegionsBuilder {

    private final int minGap;
    private final Map<String, GenomeRegions> contigMap = new LinkedHashMap<>();

    public GenomeRegionsBuilder() {
        this(1);
    }

    public GenomeRegionsBuilder(final int minGap) {
        this.minGap = minGap;
    }

    public void addPosition(@NotNull String contig, int position) {
        contigMap.computeIfAbsent(contig, x -> new GenomeRegions(contig, minGap)).addPosition(position);
    }

    public void addPosition(@NotNull GenomePosition position) {
        addPosition(position.chromosome(), position.position());
    }

    public void addRegion(@NotNull String contig, int start, int end) {
        contigMap.computeIfAbsent(contig, x -> new GenomeRegions(contig, minGap)).addRegion(start, end);
    }

    public void addRegion(@NotNull GenomeRegion region) {
        addRegion(region.chromosome(), region.start(), region.end());
    }

    @NotNull
    public List<GenomeRegion> build() {
        List<GenomeRegion> result = Lists.newArrayList();
        for (String contig : contigMap.keySet()) {
            result.addAll(contigMap.get(contig).build());
        }

        Collections.sort(result);
        return result;
    }

}

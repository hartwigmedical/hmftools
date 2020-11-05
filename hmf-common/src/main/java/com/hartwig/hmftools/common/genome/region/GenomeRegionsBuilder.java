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

    public void addPosition(GenomePosition position) {
        contigMap.computeIfAbsent(position.chromosome(), contig -> new GenomeRegions(contig, minGap)).addPosition(position.position());
    }

    public void addRegion(GenomeRegion position) {
        contigMap.computeIfAbsent(position.chromosome(), contig -> new GenomeRegions(contig, minGap))
                .addRegion(position.start(), position.end());
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

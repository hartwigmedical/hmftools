package com.hartwig.hmftools.purple.region;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

class GCAccumulator implements Consumer<GCProfile> {

    private final GenomeRegion region;
    private double totalContent;
    private int count;

    GCAccumulator(final GenomeRegion region) {
        this.region = region;
    }

    double averageGCContent() {
        return count == 0 ? 0 : totalContent / count;
    }

    @Override
    public void accept(final GCProfile gcProfile) {
        if (gcProfile.isMappable() && gcProfile.start() >= region.start() && gcProfile.end() <= region.end()) {
            count++;
            totalContent += gcProfile.gcContent();
        }
    }
}

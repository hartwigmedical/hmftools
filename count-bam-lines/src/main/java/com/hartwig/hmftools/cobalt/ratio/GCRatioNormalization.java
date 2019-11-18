package com.hartwig.hmftools.cobalt.ratio;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCount;
import com.hartwig.hmftools.common.genome.gc.GCMedianReadCountBuilder;
import com.hartwig.hmftools.common.genome.gc.GCProfile;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

class GCRatioNormalization {

    private final GCMedianReadCountBuilder medianReadCountBuilder = new GCMedianReadCountBuilder();
    private final Multimap<Chromosome, ReadCountWithGCContent> entries = ArrayListMultimap.create();

    void addPosition(@NotNull final Chromosome chromosome, @NotNull final GCProfile gcProfile, final int readCount) {
        final ReadCountWithGCContent readCountWithGCContent = new ReadCountWithGCContent(readCount, gcProfile);
        entries.put(chromosome, readCountWithGCContent);

        // TODO: TEST With/without isMappable
        if (chromosome.isAutosome() && readCountWithGCContent.isMappable() && readCount > 0) {
            medianReadCountBuilder.add(gcProfile, readCount);
        }
    }

    GCMedianReadCount gcMedianReadCount() {
        return medianReadCountBuilder.build();
    }

    @NotNull
    ListMultimap<Chromosome, ReadRatio> build(@NotNull final GCMedianReadCount gcMedianReadCount) {
        final ListMultimap<Chromosome, ReadRatio> result = ArrayListMultimap.create();
        for (Chromosome chromosome : entries.keySet()) {
            final List<ReadRatio> normalisedRatio = entries.get(chromosome).stream().map(x -> create(gcMedianReadCount, x)).collect(Collectors.toList());
            result.replaceValues(chromosome, normalisedRatio);
        }

        return result;
    }

    @NotNull
    private static ReadRatio create(@NotNull final GCMedianReadCount medians, @NotNull final ReadCountWithGCContent readCount) {
        int gcMedianCount = medians.medianReadCount(readCount.gcProfile());
        final double ratio;

        double medianNormalisation = 1.0 * medians.medianReadCount() / medians.meanReadCount();

        if (gcMedianCount == -1 || !readCount.isMappable() || gcMedianCount == 0) {
            ratio = -1;
        } else {
            ratio = medianNormalisation * readCount.readCount() / gcMedianCount;
        }

        return ImmutableReadRatio.builder().from(readCount).ratio(ratio).build();
    }

    private static class ReadCountWithGCContent implements GenomePosition {

        private final GCProfile gcProfile;
        private final int readCount;

        private ReadCountWithGCContent(final int readCount, @NotNull final GCProfile gcProfile) {
            this.readCount = readCount;
            this.gcProfile = gcProfile;
        }

        @NotNull
        @Override
        public String chromosome() {
            return gcProfile.chromosome();
        }

        @Override
        public long position() {
            return gcProfile.start();
        }

        private int readCount() {
            return readCount;
        }

        GCProfile gcProfile() {
            return gcProfile;
        }

        private boolean isMappable() {
            return gcProfile.isMappable();
        }
    }
}

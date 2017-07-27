package com.hartwig.hmftools.common.purple.ratio;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

public class NormalizedRatiosBuilder {

    private static final double MIN_MAPPABLE_PERCENTAGE = 0.85;

    private final GCMedianFactory gcMedian = new GCMedianFactory();
    private final Multimap<String, ReadCountWithGCContent> entries = ArrayListMultimap.create();

    public void addPosition(@NotNull final GCContent gcContent, @NotNull final ReadCount readCount) {
        if (gcContent.compareTo(readCount) != 0) {
            throw new IllegalArgumentException();
        }

        final ReadCountWithGCContent readCountWithGCContent = new ReadCountWithGCContent(readCount, gcContent);
        entries.put(gcContent.chromosome(), readCountWithGCContent);

        // TODO: TEST With/without ismappable
        if (!readCount.chromosome().equals("X") && !readCount.chromosome().equals("Y") && readCountWithGCContent.isMappable()) {
            gcMedian.addReadCount(readCountWithGCContent.gcContent(), readCountWithGCContent.readCount());
        }
    }

    public NormalizedRatios build() {

        final Map<Integer, GCMedians> medianCountPerGCBucket = gcMedian.medianCountPerGCBucket();
        final List<GCMedians> medians = Lists.newArrayList(medianCountPerGCBucket.values());
        Collections.sort(medians);
        final ImmutableNormalizedRatios.Builder builder = ImmutableNormalizedRatios.builder().addAllMedianReadCount(medians);

        for (String chromosome : entries.keySet()) {
            final List<ReadRatio> normalisedRatio =
                    entries.get(chromosome).stream().map(x -> create(medianCountPerGCBucket, x)).collect(Collectors.toList());
            builder.putAllNormalisedRatios(chromosome, normalisedRatio);
        }

        return builder.build();
    }

    private ReadRatio create(Map<Integer, GCMedians> medianCountPerGCBucket, ReadCountWithGCContent readCount) {
        GCMedians gcMedian = medianCountPerGCBucket.get(readCount.gcContent());
        int gcMedianCount = gcMedian == null ? -1 : gcMedian.medianCount();

        final double ratio;

        if (gcMedianCount == -1 || !readCount.isMappable() || gcMedianCount == 0) {
            ratio = -1;
        } else {
            ratio = 1.0 * readCount.readCount() / gcMedianCount;
        }

        return ImmutableReadRatio.builder().from(readCount).ratio(ratio).build();

    }

    private class ReadCountWithGCContent implements GenomePosition {

        private final GCContent gcContent;
        private final ReadCount readCount;

        private ReadCountWithGCContent(@NotNull final ReadCount readCount, @NotNull final GCContent gcContent) {
            this.readCount = readCount;
            this.gcContent = gcContent;
        }

        @NotNull
        @Override
        public String chromosome() {
            return readCount.chromosome();
        }

        @Override
        public long position() {
            return readCount.position();
        }

        private int readCount() {
            return readCount.readCount();
        }

        private int gcContent() {
            return (int) Math.round(gcContent.gcContent() * 100);
        }

        private boolean isMappable() {
            return Doubles.greaterOrEqual(gcContent.mappablePercentage(), MIN_MAPPABLE_PERCENTAGE);
        }
    }
}

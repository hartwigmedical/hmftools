package com.hartwig.hmftools.fastqstats;

import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.stream.Collectors;

import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;

import org.jetbrains.annotations.NotNull;

class FastqTracker {

    @NotNull
    private final Map<String, Map<String, FastqData>> samples;

    FastqTracker() {
        samples = ImmutableSortedMap.of();
    }

    private FastqTracker(@NotNull final Map<String, Map<String, FastqData>> samples) {
        this.samples = samples;
    }

    @NotNull
    FastqTracker addToSample(@NotNull final String sample, @NotNull final String lane, @NotNull final FastqData data) {
        final Map<String, Map<String, FastqData>> newSamples;
        if (!samples.containsKey(sample)) {
            newSamples = addToMap(samples, sample,
                    new ImmutableSortedMap.Builder<String, FastqData>(Ordering.natural()).put(lane, data).build());
        } else if (!samples.get(sample).containsKey(lane)) {
            newSamples = addToMap(samples, sample, addToMap(samples.get(sample), lane, data));
        } else {
            final FastqData currentValue = samples.get(sample).get(lane);
            final FastqData newValue = currentValue.add(data);
            newSamples = addToMap(samples, sample, addToMap(samples.get(sample), lane, newValue));
        }
        return new FastqTracker(newSamples);
    }

    @NotNull
    FastqData flowcell() {
        return samples.keySet().stream().map(this::sample).reduce(new FastqData(0, 0), FastqData::add);
    }

    @NotNull
    FastqData lane(@NotNull final String lane) {
        return samples.values().stream().map(lanes -> lanes.get(lane)).filter(Objects::nonNull).reduce(
                new FastqData(0, 0), FastqData::add);
    }

    @NotNull
    FastqData sample(@NotNull final String sample) {
        return samples.get(sample).values().stream().reduce(new FastqData(0, 0), FastqData::add);
    }

    @NotNull
    Map<String, FastqData> lanes() {
        return samples.values().stream().flatMap(sampleLanes -> sampleLanes.entrySet().stream()).collect(
                Collectors.groupingBy(Map.Entry::getKey, TreeMap::new, Collectors.mapping(Map.Entry::getValue,
                        Collectors.reducing(new FastqData(0, 0), FastqData::add))));
    }

    @NotNull
    Map<String, Map<String, FastqData>> samples() {
        return samples;
    }

    @NotNull
    private static <T> Map<String, T> addToMap(@NotNull final Map<String, T> map, @NotNull final String key,
            @NotNull final T value) {
        final Map<String, T> copy = Maps.newHashMap(map);
        copy.put(key, value);
        return ImmutableSortedMap.copyOf(copy);
    }
}

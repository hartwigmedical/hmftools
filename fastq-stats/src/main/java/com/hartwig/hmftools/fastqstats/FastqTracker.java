package com.hartwig.hmftools.fastqstats;

import java.util.HashMap;
import java.util.Map;
import java.util.Objects;
import java.util.TreeMap;
import java.util.stream.Collectors;

import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Ordering;

import org.jetbrains.annotations.NotNull;

class FastqTracker {
    private final FastqData flowcell;
    private final Map<String, Map<String, FastqData>> samples;

    FastqTracker() {
        flowcell = new FastqData(0, 0);
        samples = ImmutableSortedMap.of();
    }

    private FastqTracker(@NotNull FastqData flowcell, @NotNull Map<String, Map<String, FastqData>> samples) {
        this.flowcell = flowcell;
        this.samples = samples;
    }

    @NotNull
    FastqTracker addToFlowcell(@NotNull FastqData data) {
        return new FastqTracker(flowcell.add(data), samples);
    }

    @NotNull
    FastqTracker addToSample(@NotNull String sample, @NotNull String lane, @NotNull FastqData data) {
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
        return new FastqTracker(flowcell, newSamples);
    }

    @NotNull
    FastqData flowcell() {
        return flowcell;
    }

    @NotNull
    FastqData lane(@NotNull String lane) {
        return samples.values().stream().map(lanes -> lanes.get(lane)).filter(Objects::nonNull).reduce(
                new FastqData(0, 0), FastqData::add);
    }

    @NotNull
    FastqData sample(@NotNull String sample) {
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
    private <T> Map<String, T> addToMap(@NotNull Map<String, T> map, @NotNull String key, @NotNull T value) {
        final Map<String, T> copy = new HashMap<>(map);
        copy.put(key, value);
        return ImmutableSortedMap.copyOf(copy);
    }
}

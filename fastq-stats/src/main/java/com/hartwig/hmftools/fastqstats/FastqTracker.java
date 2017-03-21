package com.hartwig.hmftools.fastqstats;

import java.util.HashMap;
import java.util.Map;

import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.Ordering;

import org.jetbrains.annotations.NotNull;

class FastqTracker {
    private final FastqData flowcell;
    private final Map<String, FastqData> lanes;
    private final Map<String, FastqData> samples;
    private final FastqData undetermined;

    FastqTracker() {
        flowcell = new FastqData(0, 0);
        undetermined = new FastqData(0, 0);
        lanes = ImmutableSortedMap.of();
        samples = ImmutableSortedMap.of();
    }

    private FastqTracker(@NotNull FastqData flowcell, @NotNull FastqData undetermined,
            @NotNull Map<String, FastqData> lanes, @NotNull Map<String, FastqData> samples) {
        this.flowcell = flowcell;
        this.undetermined = undetermined;
        this.lanes = lanes;
        this.samples = samples;
    }

    @NotNull
    FastqTracker addToFlowcell(@NotNull FastqData data) {
        return new FastqTracker(flowcell.add(data), undetermined, lanes, samples);
    }

    @NotNull
    FastqTracker addToUndetermined(@NotNull FastqData data) {
        return new FastqTracker(flowcell, undetermined.add(data), lanes, samples);
    }

    @NotNull
    FastqTracker addToLane(@NotNull String lane, @NotNull FastqData data) {
        final Map<String, FastqData> newLanes = addToMap(lanes, lane, data);
        return new FastqTracker(flowcell, undetermined, newLanes, samples);
    }

    @NotNull
    FastqTracker addToSample(@NotNull String sample, @NotNull FastqData data) {
        final Map<String, FastqData> newSamples = addToMap(samples, sample, data);
        return new FastqTracker(flowcell, undetermined, lanes, newSamples);
    }

    @NotNull
    FastqData flowcell() {
        return flowcell;
    }

    @NotNull
    FastqData undetermined() {
        return undetermined;
    }

    @NotNull
    FastqData lane(@NotNull String lane) {
        return lanes.get(lane);
    }

    @NotNull
    FastqData sample(@NotNull String sample) {
        return samples.get(sample);
    }

    @NotNull
    Map<String, FastqData> lanes() {
        return lanes;
    }

    @NotNull
    Map<String, FastqData> samples() {
        return samples;
    }

    @NotNull
    private Map<String, FastqData> addToMap(@NotNull Map<String, FastqData> map, @NotNull String key,
            @NotNull FastqData value) {
        if (!map.containsKey(key)) {
            return new ImmutableSortedMap.Builder<String, FastqData>(Ordering.natural()).putAll(map).put(key,
                    value).build();
        } else {
            final Map<String, FastqData> copy = new HashMap<>(map);
            copy.put(key, map.get(key).add(value));
            return ImmutableSortedMap.copyOf(copy);
        }
    }
}

package com.hartwig.hmftools.fastqstats;

import java.util.Map;
import java.util.TreeMap;

import org.jetbrains.annotations.NotNull;

public class FastqTracker {
    private FastqData flowcell;
    private TreeMap<String, FastqData> lanes;
    private TreeMap<String, FastqData> samples;
    private FastqData undetermined;

    public FastqTracker() {
        flowcell = new FastqData(0, 0);
        undetermined = new FastqData(0, 0);
        lanes = new TreeMap<>();
        samples = new TreeMap<>();
    }

    public void addToFlowcell(@NotNull FastqData data) {
        flowcell = flowcell.add(data);
    }

    public void addToUndetermined(@NotNull FastqData data) {
        undetermined = undetermined.add(data);
    }

    public void addToLane(@NotNull String lane, @NotNull FastqData data) {
        if (!lanes.containsKey(lane)) {
            lanes.put(lane, data);
        } else {
            FastqData current = lanes.get(lane);
            lanes.put(lane, current.add(data));
        }
    }

    public void addToSample(@NotNull String sample, @NotNull FastqData data) {
        if (!samples.containsKey(sample)) {
            samples.put(sample, data);
        } else {
            FastqData current = samples.get(sample);
            samples.put(sample, current.add(data));
        }
    }

    @NotNull
    public FastqData getFlowcellData() {
        return flowcell;
    }

    @NotNull
    public FastqData getUndeterminedData() {
        return undetermined;
    }

    @NotNull
    public FastqData getLaneData(@NotNull String lane) {
        return lanes.get(lane);
    }

    @NotNull
    public FastqData getSampleData(@NotNull String sample) {
        return samples.get(sample);
    }

    @NotNull
    public Map<String, FastqData> getLanes() {
        return lanes;
    }

    @NotNull
    public Map<String, FastqData> getSamples() {
        return samples;
    }
}

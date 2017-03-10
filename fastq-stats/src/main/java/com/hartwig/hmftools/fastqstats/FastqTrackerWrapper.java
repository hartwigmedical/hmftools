package com.hartwig.hmftools.fastqstats;

import org.jetbrains.annotations.NotNull;

public class FastqTrackerWrapper {
    private FastqTracker tracker;

    public FastqTracker getTracker() {
        return tracker;
    }

    public FastqTrackerWrapper() {
        tracker = new FastqTracker();
    }

    public synchronized void addDataFromSampleFile(@NotNull String sampleName, @NotNull String lane,
            @NotNull FastqData data) {
        tracker = tracker.addToFlowcell(data).addToLane(lane, data).addToSample(sampleName, data);
    }

    public synchronized void addDataFromUndeterminedFile(@NotNull String lane, @NotNull FastqData data) {
        tracker = tracker.addToFlowcell(data).addToLane(lane, data).addToUndetermined(data);
    }
}

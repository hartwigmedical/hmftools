package com.hartwig.hmftools.fastqstats;

import org.jetbrains.annotations.NotNull;

class FastqTrackerWrapper {
    private FastqTracker tracker;

    FastqTrackerWrapper() {
        tracker = new FastqTracker();
    }

    synchronized FastqTracker tracker() {
        return tracker;
    }

    synchronized void addDataFromSampleFile(@NotNull String sampleName, @NotNull String lane,
            @NotNull FastqData data) {
        tracker = tracker.addToFlowcell(data).addToSample(sampleName, lane, data);
    }
}

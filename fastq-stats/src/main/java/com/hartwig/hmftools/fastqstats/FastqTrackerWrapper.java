package com.hartwig.hmftools.fastqstats;

import org.jetbrains.annotations.NotNull;

class FastqTrackerWrapper {

    @NotNull
    private FastqTracker tracker;

    FastqTrackerWrapper() {
        tracker = new FastqTracker();
    }

    @NotNull
    synchronized FastqTracker tracker() {
        return tracker;
    }

    synchronized void addDataFromSampleFile(@NotNull final String sampleName, @NotNull final String lane,
            @NotNull final FastqData data) {
        tracker = tracker.addToSample(sampleName, lane, data);
    }
}

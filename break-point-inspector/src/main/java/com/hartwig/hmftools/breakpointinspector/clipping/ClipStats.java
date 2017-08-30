package com.hartwig.hmftools.breakpointinspector.clipping;

import com.hartwig.hmftools.breakpointinspector.Location;

import org.jetbrains.annotations.NotNull;

public class ClipStats implements Comparable<ClipStats> {
    public final Location Alignment;
    public String LongestClipSequence = "";
    int Support = 1;
    public boolean Left = false;
    boolean Right = false;

    ClipStats(final Location alignment, final String sequence, final boolean left) {
        Alignment = alignment;
        LongestClipSequence = sequence;
        if (left) {
            Left = true;
        } else {
            Right = true;
        }
    }

    @Override
    public String toString() {
        return LongestClipSequence;
    }

    @Override
    public int compareTo(@NotNull final ClipStats o) {
        return LongestClipSequence.compareTo(o.LongestClipSequence);
    }
}

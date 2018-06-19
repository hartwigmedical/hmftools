package com.hartwig.hmftools.breakpointinspector.clipping;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.breakpointinspector.Location;

import org.jetbrains.annotations.NotNull;

public class ClipStats implements Comparable<ClipStats> {

    @NotNull
    public final Location alignment;
    @NotNull
    public String longestClipSequence;
    @NotNull
    public Set<String> supportingReads = Sets.newHashSet();
    public boolean left;

    ClipStats(@NotNull final Location alignment, @NotNull final String sequence, final boolean left) {
        this.alignment = alignment;
        this.longestClipSequence = sequence;
        this.left = left;
    }

    @Override
    public String toString() {
        return longestClipSequence;
    }

    @Override
    public int compareTo(@NotNull final ClipStats o) {
        return longestClipSequence.compareTo(o.longestClipSequence);
    }
}

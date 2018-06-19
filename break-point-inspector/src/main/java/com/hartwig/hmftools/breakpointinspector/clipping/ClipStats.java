package com.hartwig.hmftools.breakpointinspector.clipping;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.breakpointinspector.Location;

import org.jetbrains.annotations.NotNull;

public class ClipStats implements Comparable<ClipStats> {

    @NotNull
    private final Location alignment;
    @NotNull
    private String longestClipSequence;
    @NotNull
    private final Set<String> supportingReads = Sets.newHashSet();
    private final boolean left;

    ClipStats(@NotNull final Location alignment, @NotNull final String sequence, final boolean left) {
        this.alignment = alignment;
        this.longestClipSequence = sequence;
        this.left = left;
    }

    @NotNull
    public Location alignment() {
        return alignment;
    }

    void newLongestClipSequence(@NotNull final String sequence) {
        this.longestClipSequence = sequence;
    }

    @NotNull
    public String longestClipSequence() {
        return longestClipSequence;
    }

    void addSupportingRead(@NotNull String supportingRead) {
        supportingReads.add(supportingRead);
    }

    @NotNull
    public Set<String> supportingReads() {
        return supportingReads;
    }

    public boolean left() {
        return left;
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

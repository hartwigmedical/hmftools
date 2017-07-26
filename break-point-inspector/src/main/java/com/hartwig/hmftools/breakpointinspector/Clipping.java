package com.hartwig.hmftools.breakpointinspector;

import java.util.Comparator;
import java.util.List;
import java.util.SortedSet;
import java.util.stream.Collectors;

import com.google.common.collect.TreeMultimap;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

class ClipInfo {
    Location Alignment;
    int Length = 0;
    String Sequence = "";
    boolean HardClipped = false;
    boolean Left = false;
    boolean Right = false;
}

class ClipStats implements Comparable<ClipStats> {
    final Location Alignment;
    String LongestClipSequence = "";
    int Support = 1;

    ClipStats(final Location alignment, final String sequence) {
        Alignment = alignment;
        LongestClipSequence = sequence;
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

class Clipping {

    private TreeMultimap<Location, ClipStats> LocationMap = TreeMultimap.create(Comparator.naturalOrder(), Comparator.naturalOrder());

    void add(@Nullable final ClipInfo clip) {

        if (clip == null) {
            return;
        }

        boolean found = false;
        final SortedSet<ClipStats> existing = LocationMap.removeAll(clip.Alignment);
        for (final ClipStats c : existing) {
            if (clip.Left) {
                if (c.LongestClipSequence.length() > clip.Sequence.length()) {
                    if (c.LongestClipSequence.endsWith(clip.Sequence)) {
                        c.Support++;
                        found = true;
                    }
                } else if (clip.Sequence.endsWith(c.LongestClipSequence)) {
                    c.LongestClipSequence = clip.Sequence;
                    c.Support++;
                    found = true;
                }
            } else if (clip.Right) {
                if (c.LongestClipSequence.length() > clip.Sequence.length()) {
                    if (c.LongestClipSequence.startsWith(clip.Sequence)) {
                        c.Support++;
                        found = true;
                    }
                } else if (clip.Sequence.startsWith(c.LongestClipSequence)) {
                    c.LongestClipSequence = clip.Sequence;
                    c.Support++;
                    found = true;
                }
            }

            LocationMap.put(clip.Alignment, c);
        }

        if (!found) {
            LocationMap.put(clip.Alignment, new ClipStats(clip.Alignment, clip.Sequence));
        }
    }

    List<ClipStats> getSequences() {
        return LocationMap.values().stream().sorted((a, b) -> Integer.compare(b.Support, a.Support)).collect(Collectors.toList());
    }

    static ClipInfo getLeftClip(final SAMRecord read) {
        final ClipInfo result = new ClipInfo();
        result.Left = true;
        final Cigar cigar = read.getCigar();
        if (cigar.isEmpty()) {
            return null;
        }
        switch (cigar.getFirstCigarElement().getOperator()) {
            case H:
                result.Alignment = Location.fromSAMRecord(read, true).add(-1);
                result.Length = cigar.getFirstCigarElement().getLength();
                result.HardClipped = true;
                return result;
            case S:
                result.Alignment = Location.fromSAMRecord(read, true).add(-1);
                result.Length = cigar.getFirstCigarElement().getLength();
                result.Sequence = read.getReadString().substring(0, result.Length);
                return result;
        }
        return null;
    }

    static ClipInfo getRightClip(final SAMRecord read) {
        final ClipInfo result = new ClipInfo();
        result.Right = true;
        final Cigar cigar = read.getCigar();
        if (cigar.isEmpty()) {
            return null;
        }
        switch (cigar.getLastCigarElement().getOperator()) {
            case H:
                result.Alignment = Location.fromSAMRecord(read, false).add(1);
                result.Length = cigar.getLastCigarElement().getLength();
                result.HardClipped = true;
                return result;
            case S:
                result.Alignment = Location.fromSAMRecord(read, false).add(1);
                result.Length = cigar.getLastCigarElement().getLength();
                result.Sequence = read.getReadString().substring(read.getReadLength() - result.Length);
                return result;
        }
        return null;
    }
}

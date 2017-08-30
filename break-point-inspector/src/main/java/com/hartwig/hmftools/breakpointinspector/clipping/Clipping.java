package com.hartwig.hmftools.breakpointinspector.clipping;

import java.util.Comparator;
import java.util.List;
import java.util.SortedSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.breakpointinspector.Location;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

public class Clipping {

    private TreeMultimap<Location, ClipStats> LocationMap = TreeMultimap.create(Comparator.naturalOrder(), Comparator.naturalOrder());

    public void add(@Nullable final ClipInfo clip) {

        if (clip == null) {
            return;
        }

        boolean found = false;
        final SortedSet<ClipStats> existing = LocationMap.removeAll(clip.Alignment);
        for (final ClipStats c : existing) {
            if (clip.Left != c.Left) {
                continue;
            }
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
            LocationMap.put(clip.Alignment, new ClipStats(clip.Alignment, clip.Sequence, clip.Left));
        }
    }

    public List<ClipStats> getSequences() {
        return LocationMap.values().stream().sorted((a, b) -> Integer.compare(b.Support, a.Support)).collect(Collectors.toList());
    }

    public List<ClipStats> getSequencesAt(final Location location) {
        return Lists.newArrayList(LocationMap.get(location));
    }

    public static ClipInfo getLeftClip(final SAMRecord read) {
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

    public static ClipInfo getRightClip(final SAMRecord read) {
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

    public static Stream<ClipInfo> getClips(final SAMRecord read) {
        return Stream.of(getLeftClip(read), getRightClip(read));
    }
}

package com.hartwig.hmftools.breakpointinspector.clipping;

import java.util.Comparator;
import java.util.List;
import java.util.SortedSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.breakpointinspector.Location;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

public class Clipping {

    @NotNull
    private final TreeMultimap<Location, ClipStats> locationMap = TreeMultimap.create(Comparator.naturalOrder(), Comparator.naturalOrder());

    public void add(@Nullable final ClipInfo clipInfo) {
        if (clipInfo == null) {
            return;
        }

        boolean found = false;
        final SortedSet<ClipStats> existing = locationMap.removeAll(clipInfo.alignment());
        for (final ClipStats clipStats : existing) {
            if (clipInfo.left() != clipStats.left) {
                continue;
            }

            boolean addToExisting = false;
            if (clipInfo.left()) {
                if (clipStats.longestClipSequence.length() > clipInfo.sequence().length()) {
                    if (clipStats.longestClipSequence.endsWith(clipInfo.sequence())) {
                        addToExisting = true;
                    }
                } else if (clipInfo.sequence().endsWith(clipStats.longestClipSequence)) {
                    clipStats.longestClipSequence = clipInfo.sequence();
                    addToExisting = true;
                }
            } else if (clipInfo.right()) {
                if (clipStats.longestClipSequence.length() > clipInfo.sequence().length()) {
                    if (clipStats.longestClipSequence.startsWith(clipInfo.sequence())) {
                        addToExisting = true;
                    }
                } else if (clipInfo.sequence().startsWith(clipStats.longestClipSequence)) {
                    clipStats.longestClipSequence = clipInfo.sequence();
                    addToExisting = true;
                }
            }

            if (addToExisting) {
                clipStats.supportingReads.add(clipInfo.record().getReadName());
                found = true;
            }
            locationMap.put(clipInfo.alignment(), clipStats);
        }

        if (!found) {
            final ClipStats stats = new ClipStats(clipInfo.alignment(), clipInfo.sequence(), clipInfo.left());
            stats.supportingReads.add(clipInfo.record().getReadName());
            locationMap.put(clipInfo.alignment(), stats);
        }
    }

    @NotNull
    public List<ClipStats> getSequences() {
        return locationMap.values()
                .stream()
                .sorted((a, b) -> Integer.compare(b.supportingReads.size(), a.supportingReads.size()))
                .collect(Collectors.toList());
    }

    @NotNull
    public List<ClipStats> sequencesAt(@NotNull final Location location) {
        return Lists.newArrayList(locationMap.get(location));
    }

    @Nullable
    public static ClipInfo leftClip(@NotNull final SAMRecord read) {
        final Cigar cigar = read.getCigar();
        if (cigar.isEmpty()) {
            return null;
        }
        Location alignment;
        int length;
        boolean hardClipped = false;
        String sequence = "";
        switch (cigar.getFirstCigarElement().getOperator()) {
            case H:
                alignment = Location.fromSAMRecord(read, true).add(-1);
                length = cigar.getFirstCigarElement().getLength();
                hardClipped = true;
                break;
            case S:
                alignment = Location.fromSAMRecord(read, true).add(-1);
                length = cigar.getFirstCigarElement().getLength();
                sequence = read.getReadString().substring(0, length);
                break;
            default:
                return null;
        }

        return ImmutableClipInfo.of(read, alignment, length, sequence, hardClipped, true, false);
    }

    @Nullable
    public static ClipInfo rightClip(@NotNull final SAMRecord read) {
        final Cigar cigar = read.getCigar();
        if (cigar.isEmpty()) {
            return null;
        }
        Location alignment;
        int length;
        boolean hardClipped = false;
        String sequence = "";
        switch (cigar.getLastCigarElement().getOperator()) {
            case H:
                alignment = Location.fromSAMRecord(read, false).add(1);
                length = cigar.getLastCigarElement().getLength();
                hardClipped = true;
                break;
            case S:
                alignment = Location.fromSAMRecord(read, false).add(1);
                length = cigar.getLastCigarElement().getLength();
                sequence = read.getReadString().substring(read.getReadLength() - length);
                break;
            default:
                return null;
        }

        return ImmutableClipInfo.of(read, alignment, length, sequence, hardClipped, false, true);
    }

    @NotNull
    public static Stream<ClipInfo> clips(@NotNull final SAMRecord read) {
        return Stream.of(leftClip(read), rightClip(read));
    }
}

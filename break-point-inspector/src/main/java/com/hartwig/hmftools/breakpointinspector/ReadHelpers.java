package com.hartwig.hmftools.breakpointinspector;

import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;

import static com.hartwig.hmftools.breakpointinspector.Util.*;

class ReadHelpers {

    static boolean readIntersectsLocation(final SAMRecord read, final Location location) {
        assert !read.getReadUnmappedFlag();
        return read.getReferenceIndex() == location.ReferenceIndex && read.getAlignmentStart() <= location.Position
                && read.getAlignmentEnd() >= location.Position;
    }

    static boolean pairStraddlesLocation(final SAMRecord read, final Location location) {
        assert !read.getReadUnmappedFlag() && !read.getMateUnmappedFlag();
        if (read.getInferredInsertSize() > 0)
            return read.getAlignmentStart() <= location.Position
                    && (read.getAlignmentStart() + read.getInferredInsertSize()) >= location.Position;
        else
            return read.getMateAlignmentStart() <= location.Position
                    && (read.getMateAlignmentStart() - read.getInferredInsertSize()) >= location.Position;
    }

    static boolean isClipped(final SAMRecord read) {
        return read.getCigar().isClipped();
    }

    private static ClipInfo getLeftClip(final SAMRecord read) {
        final ClipInfo result = new ClipInfo();
        result.Left = true;
        final Cigar cigar = read.getCigar();
        if (cigar.isEmpty())
            return null;
        switch (cigar.getFirstCigarElement().getOperator()) {
            case H:
                result.Alignment = Location.fromSAMRecord(read, true);
                result.Length = cigar.getFirstCigarElement().getLength();
                result.HardClipped = true;
                return result;
            case S:
                result.Alignment = Location.fromSAMRecord(read, true);
                result.Length = cigar.getFirstCigarElement().getLength();
                result.Sequence = read.getReadString().substring(0, result.Length);
                return result;
        }
        return null;
    }

    private static ClipInfo getRightClip(final SAMRecord read) {
        final ClipInfo result = new ClipInfo();
        result.Right = true;
        final Cigar cigar = read.getCigar();
        if (cigar.isEmpty())
            return null;
        switch (cigar.getLastCigarElement().getOperator()) {
            case H:
                result.Alignment = Location.fromSAMRecord(read, false);
                result.Length = cigar.getLastCigarElement().getLength();
                result.HardClipped = true;
                return result;
            case S:
                result.Alignment = Location.fromSAMRecord(read, false);
                result.Length = cigar.getLastCigarElement().getLength();
                result.Sequence = read.getReadString().substring(read.getReadLength() - result.Length);
                return result;
        }
        return null;
    }

    static List<ClipInfo> getClips(final SAMRecord read) {
        return Stream.of(getLeftClip(read), getRightClip(read)).filter(Objects::nonNull).collect(Collectors.toList());
    }

    static Region determineRegion(final SAMRecord read, final Location location1,
            final Location location2) {
        if (location1.ReferenceIndex != location2.ReferenceIndex) {
            if (read.getReferenceIndex() == location1.ReferenceIndex)
                return Region.BP1;
            else if (read.getReferenceIndex() == location2.ReferenceIndex)
                return Region.BP2;
        } else if (read.getReferenceIndex() == location1.ReferenceIndex) {
            final int distance1 = Math.min(Math.abs(read.getAlignmentStart() - location1.Position),
                    Math.abs(read.getAlignmentEnd() - location1.Position));
            final int distance2 = Math.min(Math.abs(read.getAlignmentStart() - location2.Position),
                    Math.abs(read.getAlignmentEnd() - location2.Position));

            if (distance1 < 500 && distance2 < 500) {
                // TODO: handle uncertainty
            } else if (distance1 < distance2)
                return Region.BP1;
            else
                return Region.BP2;
        }
        return Region.OTHER;
    }

    static boolean isMate(final SAMRecord read, final SAMRecord mate) {
        return read.getReadName().equals(mate.getReadName()) && read.getMateReferenceIndex().equals(
                mate.getReferenceIndex()) && Math.abs(read.getMateAlignmentStart() - mate.getAlignmentStart()) <= 1;
    }
}

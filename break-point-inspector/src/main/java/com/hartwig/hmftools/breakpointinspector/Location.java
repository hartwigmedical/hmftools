package com.hartwig.hmftools.breakpointinspector;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;

public class Location implements Comparable<Location> {

    @NotNull
    private final String referenceName;
    private final int referenceIndex;
    private final int position;

    public Location(@NotNull final String referenceName, final int referenceIndex, final int position) {
        this.referenceName = referenceName;
        this.referenceIndex = referenceIndex;
        this.position = position;
    }

    public int referenceIndex() {
        return referenceIndex;
    }

    public int position() {
        return position;
    }

    @NotNull
    static Location parseFromVariant(@NotNull VariantContext variant, @NotNull SAMSequenceDictionary dictionary) {
        final String variantLocation = variant.getContig() + ":" + Integer.toString(variant.getStart());
        return Location.parseLocationString(variantLocation, dictionary);
    }

    @NotNull
    static Location parseLocationString(final String location, final SAMSequenceDictionary dictionary) {
        final String[] split = location.split(":");
        if (split.length != 2) {
            throw new RuntimeException(location + " is not a valid location string");
        }

        final String chromosome = split[0];

        int position;
        try {
            position = Integer.parseInt(split[1]);
        } catch (NumberFormatException e) {
            throw new RuntimeException(location + " is not a valid location string");
        }

        // NERA: Query the position
        int referenceIndex = dictionary.getSequenceIndex(chromosome);
        if (referenceIndex < 0) {
            if (!chromosome.startsWith("chr")) {
                referenceIndex = dictionary.getSequenceIndex("chr" + chromosome);
            } else {
                referenceIndex = dictionary.getSequenceIndex(chromosome.substring(3));
            }
        }
        if (referenceIndex < 0) {
            throw new RuntimeException(chromosome + " is not in the BAM");
        }

        return new Location(chromosome, referenceIndex, position);
    }

    @NotNull
    static Location fromSAMRecord(@NotNull final SAMRecord record) {
        return fromSAMRecord(record, true);
    }

    @NotNull
    public static Location fromSAMRecord(final SAMRecord record, boolean alignmentStart) {
        return new Location(record.getReferenceName(),
                record.getReferenceIndex(),
                alignmentStart ? record.getAlignmentStart() : record.getAlignmentEnd());
    }

    @NotNull
    public Location add(int delta) {
        return new Location(referenceName, referenceIndex, position + delta);
    }

    @NotNull
    Location withNewPosition(int newPosition) {
        return new Location(referenceName, referenceIndex, newPosition);
    }

    boolean sameChromosomeAs(final Location other) {
        return other != null && other.referenceIndex == referenceIndex;
    }

    @Override
    public int hashCode() {
        return Objects.hash(referenceIndex, position);
    }

    @Override
    public boolean equals(final Object obj) {
        if (obj == null) {
            return false;
        }
        if (!(obj instanceof Location)) {
            return false;
        }

        Location other = (Location) obj;
        return referenceIndex == other.referenceIndex && position == other.position;
    }

    @Override
    public int compareTo(@NotNull final Location o) {
        final int comp1 = Integer.compare(referenceIndex, o.referenceIndex);
        if (comp1 == 0) {
            return Integer.compare(position, o.position);
        }
        return comp1;
    }

    @Override
    public String toString() {
        return referenceName + ":" + position;
    }
}

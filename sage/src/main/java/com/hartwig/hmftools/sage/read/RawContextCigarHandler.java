package com.hartwig.hmftools.sage.read;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.sam.CigarHandler;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class RawContextCigarHandler implements CigarHandler {

    private final VariantHotspot variant;
    private final boolean isInsert;
    private final boolean isDelete;
    private final boolean isSNV;

    private RawContext result;

    RawContextCigarHandler(final VariantHotspot variant) {
        this.variant = variant;
        this.isInsert = variant.ref().length() < variant.alt().length();
        this.isDelete = variant.ref().length() > variant.alt().length();
        this.isSNV = variant.ref().length() == variant.alt().length();

    }

    public RawContext result() {
        return result;
    }

    @Override
    public void handleLeftSoftClip(@NotNull final SAMRecord record, @NotNull final CigarElement element) {

        if (variant.position() < record.getAlignmentStart()) {
            int readIndex = record.getReadPositionAtReferencePosition(record.getAlignmentStart()) - 1 - record.getAlignmentStart()
                    + (int) variant.position() - variant.alt().length() + variant.ref().length();
            result = RawContext.clipped(readIndex);
        }
    }

    @Override
    public void handleRightSoftClip(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
            final int refPosition) {
        if (result != null) {
            return;
        }

        long refPositionEnd = refPosition + element.getLength() - 1;
        if (refPositionEnd < variant.position()) {
            throw new IllegalStateException("Variant is after record");
        }

        if (variant.position() >= refPosition && variant.position() <= refPositionEnd) {
            int alignmentEnd = record.getAlignmentEnd();
            int actualIndex = record.getReadPositionAtReferencePosition(alignmentEnd) - 1 - alignmentEnd + (int) variant.position();
            result = RawContext.clipped(actualIndex);
        }

    }

    @Override
    public void handleAlignment(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
            final int refPosition) {
        if (result != null) {
            return;
        }

        long refPositionEnd = refPosition + element.getLength() - 1;
        if (refPosition <= variant.position() && variant.position() <= refPositionEnd) {
            int readIndexOffset = (int) (variant.position() - refPosition);
            int variantReadIndex = readIndex + readIndexOffset;

            int baseQuality = record.getBaseQualities()[variantReadIndex];
            boolean altSupport = isSNV && refPositionEnd >= variant.end() && matchesString(record, variantReadIndex, variant.alt());
            boolean refSupport = !altSupport && matchesString(record, variantReadIndex, variant.ref().substring(0, 1));
            result = RawContext.snv(variantReadIndex, altSupport, refSupport, baseQuality);
        }

    }

    @Override
    public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPosition) {
        if (result != null) {
            return;
        }

        if (refPosition == variant.position()) {
            boolean altSupport = isInsert && e.getLength() == variant.alt().length() - 1 && matchesString(record, readIndex, variant.alt());
            int baseQuality = altSupport ? baseQuality(readIndex, record, variant.alt().length()) : 0;
            result = RawContext.indel(readIndex, altSupport, baseQuality);
        }

    }

    @Override
    public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement e, final int readIndex, final int refPosition) {
        if (result != null) {
            return;
        }

        int refPositionEnd = refPosition + e.getLength();
        if (refPosition == variant.position()) {
            boolean altSupport = isDelete && e.getLength() == variant.ref().length() - 1 && matchesString(record,
                    readIndex,
                    variant.ref().substring(0, 1));
            int baseQuality = altSupport ? baseQuality(readIndex, record, 2) : 0;
            result = RawContext.indel(readIndex, altSupport, baseQuality);
        } else if (refPositionEnd >= variant.position()) {
            result = RawContext.inDelete(readIndex);
        }

    }

    @Override
    public void handleSkippedReference(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
            final int refPosition) {
        handleDelete(record, element, readIndex, refPosition);
    }

    private static boolean matchesString(SAMRecord record, int index, String expected) {
        return new String(record.getReadBases(), index, expected.length()).equals(expected);
    }

    private int baseQuality(int readIndex, SAMRecord record, int length) {
        int maxIndex = Math.min(readIndex + length, record.getBaseQualities().length) - 1;
        int quality = Integer.MAX_VALUE;
        for (int i = readIndex; i <= maxIndex; i++) {
            quality = Math.min(quality, record.getBaseQualities()[i]);
        }
        return quality;
    }
}

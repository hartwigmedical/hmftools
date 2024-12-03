package com.hartwig.hmftools.bamtools.remapper;

import java.util.Objects;

import com.hartwig.hmftools.common.bam.SamRecordUtils;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class RecordPair
{
    public final @NotNull SAMRecord left;
    public final @NotNull SAMRecord right;

    public RecordPair(@NotNull final SAMRecord left, @NotNull final SAMRecord right)
    {
        if (!left.getReadName().equals(right.getReadName())) {
            throw new IllegalArgumentException("Left read name does not match right read name.");
        }
        if (!SamRecordUtils.firstInPair(left)) {
            throw new IllegalArgumentException("Left read is not first in pair.");
        }
        if (!SamRecordUtils.secondInPair(right)) {
            throw new IllegalArgumentException("Right read is not second in pair.");
        }
        this.left = left;
        this.right = right;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final RecordPair that = (RecordPair) o;
        return Objects.equals(left, that.left) && Objects.equals(right, that.right);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(left, right);
    }
}

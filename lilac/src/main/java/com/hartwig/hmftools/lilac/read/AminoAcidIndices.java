/*
 * Decompiled with CFR 0.151.
 *
 * Could not load the following classes:
 *  kotlin.Metadata
 *  kotlin.ranges.IntRange
 *  org.jetbrains.annotations.NotNull
 */
package com.hartwig.hmftools.lilac.read;

import kotlin.Metadata;
import kotlin.ranges.IntRange;

import org.jetbrains.annotations.NotNull;

@Metadata(mv = { 1, 1, 18 },
          bv = { 1, 0, 3 },
          k = 1,
          d1 = { "\u0000\u001a\n\u0002\u0018\u0002\n\u0002\u0010\u0000\n\u0002\b\u0002\n\u0002\u0018\u0002\n\u0000\n\u0002\u0010\b\n\u0002\b\u0002\b\u00c6\u0002\u0018\u00002\u00020\u0001B\u0007\b\u0002\u00a2\u0006\u0002\u0010\u0002J\u0016\u0010\u0003\u001a\u00020\u00042\u0006\u0010\u0005\u001a\u00020\u00062\u0006\u0010\u0007\u001a\u00020\u0006\u00a8\u0006\b" },
          d2 = { "Lcom/hartwig/hmftools/lilac/read/AminoAcidIndices;", "", "()V", "indices", "Lkotlin/ranges/IntRange;", "nucStartIndex",
                  "", "nucEndIndex", "lilac" })
public final class AminoAcidIndices
{
    public static final AminoAcidIndices INSTANCE;

    @NotNull
    public final IntRange indices(int nucStartIndex, int nucEndIndex)
    {
        int start = nucStartIndex / 3 + (nucStartIndex % 3 == 0 ? 0 : 1);
        int end = (nucEndIndex + 1) / 3 - 1;
        return new IntRange(start, end);
    }

    private AminoAcidIndices()
    {
    }

    static
    {
        AminoAcidIndices aminoAcidIndices;
        INSTANCE = aminoAcidIndices = new AminoAcidIndices();
    }
}

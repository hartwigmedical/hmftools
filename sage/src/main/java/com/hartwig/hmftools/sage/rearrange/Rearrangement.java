package com.hartwig.hmftools.sage.rearrange;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class Rearrangement
{
    private static final int MAX_MNV_LENGTH = 2;

    public static List<VariantHotspot> rearrangeInsert(int position, int indelLength, int readIndex, byte[] readBases, int refIndex,
            byte[] refBases)
    {
        final List<VariantHotspot> left = moveLeft(false, position, indelLength, readIndex, readBases, refIndex, refBases);
        if(!left.isEmpty())
        {
            return left;
        }
        return moveRight(false, position, indelLength, readIndex, readBases, refIndex, refBases);
    }

    public static List<VariantHotspot> rearrangeDelete(int position, int indelLength, int readIndex, byte[] readBases, int refIndex,
            byte[] refBases)
    {
        final List<VariantHotspot> left = moveLeft(true, position, indelLength, refIndex, refBases, readIndex, readBases);
        if(!left.isEmpty())
        {
            return left;
        }
        return moveRight(true, position, indelLength, refIndex, refBases, readIndex, readBases);
    }

    @NotNull
    public static List<VariantHotspot> moveLeft(boolean delete, int position, int indelLength, int readIndex, byte[] readBases,
            int refIndex, byte[] refBases)
    {
        int absIndelLength = Math.abs(indelLength);
        int minReadIndex = Math.max(0, readIndex - absIndelLength);
        int maxRefIndex = Math.max(0, refIndex - absIndelLength);

        int indelIndex = readIndex + absIndelLength - 1;

        int sameBases = 0;
        while(refIndex >= maxRefIndex && indelIndex >= minReadIndex && readBases[indelIndex] == refBases[refIndex])
        {
            sameBases++;
            indelIndex--;
            refIndex--;
        }

        int snvLength = 0;
        while(refIndex >= maxRefIndex && indelIndex >= minReadIndex && readBases[indelIndex] != refBases[refIndex])
        {
            snvLength++;
            indelIndex--;
            refIndex--;
        }

        int sameBases2 = 0;
        while(refIndex >= maxRefIndex && indelIndex >= minReadIndex && readBases[indelIndex] == refBases[refIndex])
        {
            sameBases2++;
            indelIndex--;
            refIndex--;
        }

        if(sameBases == 0 || snvLength == 0 || sameBases2 == 0 || snvLength > MAX_MNV_LENGTH)
        {
            return Collections.emptyList();
        }

        int indelShift = sameBases + snvLength + sameBases2;

        final List<VariantHotspot> result = Lists.newArrayList();

        final String indelRef = new String(refBases, refIndex, 1);
        final String indelAlt = new String(readBases, readIndex - indelShift, absIndelLength);

        VariantHotspot indel = ImmutableVariantHotspotImpl.builder()
                .chromosome("1")
                .position(position - indelShift)
                .ref(delete ? indelAlt : indelRef)
                .alt(delete ? indelRef : indelAlt)
                .build();
        result.add(indel);

        final String snvRef = new String(refBases, refIndex + sameBases2 + 1, snvLength);
        final String snvAlt = new String(readBases, refIndex + sameBases2 + 1 + absIndelLength - 1, snvLength);

        VariantHotspot snv = ImmutableVariantHotspotImpl.builder()
                .chromosome("1")
                .position(position - indelShift + sameBases2 + 1 + (delete ? absIndelLength - 1 : 0))
                .ref(delete ? snvAlt : snvRef)
                .alt(delete ? snvRef : snvAlt)
                .build();
        result.add(snv);

        return result;
    }

    @NotNull
    public static List<VariantHotspot> moveRight(boolean delete, int position, int indelLength, int readIndex, byte[] readBases,
            int refIndex, byte[] refBases)
    {
        int absIndelLength = Math.abs(indelLength);
        int maxReadIndex = Math.min(readBases.length - 1, readIndex + 2 * absIndelLength + 1);
        int maxRefIndex = Math.min(refBases.length - 1, refIndex + 2 * absIndelLength + 1);

        // Assume everything is pointing before indel
        position++;
        readIndex++;
        refIndex++;

        int snvPosition = position;
        int snvReadIndex = readIndex;
        int snvRefIndex = refIndex;
        while(readIndex <= maxReadIndex && refIndex <= maxRefIndex && readBases[readIndex] != refBases[refIndex])
        {
            readIndex++;
            refIndex++;
            position++;
        }

        int gapPosition = position;
        while(readIndex <= maxReadIndex && refIndex <= maxRefIndex && readBases[readIndex] == refBases[refIndex])
        {
            readIndex++;
            refIndex++;
            position++;
        }

        int snvLength = gapPosition - snvPosition;
        int gapLength = position - gapPosition;

        if(snvLength == 0 || gapLength == 0 || snvLength > MAX_MNV_LENGTH)
        {
            return Collections.emptyList();
        }

        final String snvRef = new String(refBases, snvRefIndex, snvLength);
        final String snvAlt = new String(readBases, snvReadIndex, snvLength);

        final List<VariantHotspot> result = Lists.newArrayList();
        VariantHotspot snv = ImmutableVariantHotspotImpl.builder()
                .chromosome("1")
                .position(snvPosition)
                .ref(delete ? snvAlt : snvRef)
                .alt(delete ? snvRef : snvAlt)
                .build();
        result.add(snv);

        final String indelRef = new String(refBases, snvRefIndex + snvLength - 1 + gapLength, 1);
        final String indelAlt = new String(readBases, snvReadIndex + snvLength - 1 + gapLength, indelLength);

        VariantHotspot indel = ImmutableVariantHotspotImpl.builder()
                .chromosome("1")
                .position(snvPosition + snvLength - 1 + gapLength)
                .ref(delete ? indelAlt : indelRef)
                .alt(delete ? indelRef : indelAlt)
                .build();
        result.add(indel);

        return result;
    }
}

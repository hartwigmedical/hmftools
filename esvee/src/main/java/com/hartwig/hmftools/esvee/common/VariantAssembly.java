package com.hartwig.hmftools.esvee.common;

import java.util.Objects;

import com.hartwig.hmftools.esvee.models.AlignedAssembly;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;

public class VariantAssembly
{
    public final AlignedAssembly Assembly;
    @Nullable
    public final Cigar LeftAnchorCigar;
    @Nullable
    public final Cigar RightAnchorCigar;
    public final int LeftCigarLength;
    public final int RightCigarLength;
    /**
     * Position in assembly
     */
    public final int LeftPosition;
    /**
     * Position in assembly
     */
    public final int RightPosition;
    public final int LeftOverhang;
    public final int RightOverhang;

    public VariantAssembly(final AlignedAssembly assembly,
            @Nullable final Cigar leftAnchorCigar,
            final int leftCigarLength,
            final int leftPosition,
            final int leftOverhang,
            @Nullable final Cigar rightAnchorCigar,
            final int rightCigarLength,
            final int rightPosition,
            final int rightOverhang)
    {
        Assembly = assembly;
        LeftAnchorCigar = leftAnchorCigar;
        LeftCigarLength = leftCigarLength;
        LeftPosition = leftPosition;
        LeftOverhang = leftOverhang;
        RightAnchorCigar = rightAnchorCigar;
        RightCigarLength = rightCigarLength;
        RightPosition = rightPosition;
        RightOverhang = rightOverhang;
    }

    public static VariantAssembly create(final AlignedAssembly assembly,
            @Nullable final Pair<Cigar, Integer> leftAnchor,
            final int leftPosition, final int leftOverhang,
            @Nullable final Pair<Cigar, Integer> rightAnchor,
            final int rightPosition, final int rightOverhang)
    {
        return new VariantAssembly(
                assembly,
                leftAnchor == null ? null : leftAnchor.getKey(),
                leftAnchor == null ? 0 : leftAnchor.getValue(),
                leftPosition, leftOverhang,
                rightAnchor == null ? null : rightAnchor.getKey(),
                rightAnchor == null ? 0 : rightAnchor.getValue(),
                rightPosition, rightOverhang
        );
    }

    public VariantAssembly reverse()
    {
        return new VariantAssembly(Assembly,
                RightAnchorCigar,
                RightCigarLength,
                RightPosition,
                RightOverhang,
                LeftAnchorCigar,
                LeftCigarLength,
                LeftPosition,
                LeftOverhang);
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
        {
            return true;
        }
        else if(o == null || getClass() != o.getClass())
        {
            return false;
        }

        final VariantAssembly that = (VariantAssembly) o;
        return LeftCigarLength == that.LeftCigarLength && LeftPosition == that.LeftPosition
                && RightCigarLength == that.RightCigarLength && RightPosition == that.RightPosition
                && Assembly.equals(that.Assembly)
                && Objects.equals(LeftAnchorCigar, that.LeftAnchorCigar) && Objects.equals(RightAnchorCigar, that.RightAnchorCigar);
    }

    @Override
    public int hashCode()
    {
        int result = Assembly.hashCode();
        result = 31 * result + (LeftAnchorCigar != null ? LeftAnchorCigar.hashCode() : 0);
        result = 31 * result + (RightAnchorCigar != null ? RightAnchorCigar.hashCode() : 0);
        result = 31 * result + LeftCigarLength;
        result = 31 * result + RightCigarLength;
        result = 31 * result + LeftPosition;
        result = 31 * result + RightPosition;
        return result;
    }
}

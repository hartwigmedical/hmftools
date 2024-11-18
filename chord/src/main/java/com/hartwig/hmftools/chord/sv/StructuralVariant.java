package com.hartwig.hmftools.chord.sv;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class StructuralVariant
{
    public final VariantContext Context;
    public final VariantAltInsertCoords ParsedAltString;

    public final String Id;

    public final String RefChromosome;
    public final int RefPosition;
    public final Orientation RefOrientation;

    public String AltChromosome;
    public int AltPosition;
    public Orientation AltOrientation;

    public StructuralVariantType Type;
    @Nullable public Integer Length;

    private StructuralVariant(VariantContext context, String id, String refChromosome, int refPosition, String refSequence, String altString)
    {
        Context = context;
        Id = id;

        ParsedAltString = VariantAltInsertCoords.fromRefAlt(altString, refSequence);

        RefChromosome = refChromosome;
        RefPosition = refPosition;
        RefOrientation = ParsedAltString.Orient;

        AltChromosome = (!ParsedAltString.OtherChromsome.isEmpty()) ?
                ParsedAltString.OtherChromsome :
                null;
        AltPosition = ParsedAltString.OtherPosition;
        AltOrientation = ParsedAltString.OtherOrient;

        Type = getType();
        Length = getLength();
    }

    public StructuralVariant(VariantContext context)
    {
        this(
                context,
                context.getID(),
                context.getContig(),
                context.getStart(),
                context.getReference().getBaseString(),
                context.getAlternateAllele(0).getDisplayString()
        );
    }

    @VisibleForTesting
    public StructuralVariant(String refChromosome, int refPosition, String refSequence, String altString)
    {
        this(null, null, refChromosome, refPosition, refSequence, altString);
    }

    private StructuralVariantType getType()
    {
        if(AltChromosome == null)
            return SGL;

        if(!RefChromosome.equals(AltChromosome))
            return BND;

        int posDiff = Math.abs(AltPosition - RefPosition);
        // if(posDiff == 1 && hasInsertedBases) // This is the definition ESVEE uses for insertions
        if(posDiff == 1) // However, CHORD was training using this definition
            return INS;

        if(RefOrientation == AltOrientation)
            return INV;

        boolean refPosIsLower = RefPosition < AltPosition;

        if(refPosIsLower)
            return RefOrientation.isForward() ? DEL : DUP;
        else
            return RefOrientation.isReverse() ? DEL : DUP;
    }

    private @Nullable Integer getLength()
    {
        if(Type==DEL || Type==DUP || Type == INV || Type == INS)
            return Math.abs(AltPosition - RefPosition);
        else
            return null;
    }

    public String toString()
    {
        return String.format("id(%s) coords(%s:%d:%d)", Id, Context.getContig(), RefPosition, RefOrientation.asByte());
    }

}

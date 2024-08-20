package com.hartwig.hmftools.esvee.utils.vcfcompare;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.HOMSEQ;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SVTYPE;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.EVENT_TYPE;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.esvee.common.CommonUtils.formSvType;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;

import htsjdk.variant.variantcontext.VariantContext;

class VariantBreakend
{
    public final VariantContext Context;
    public final VariantAltInsertCoords AltCoords;

    public final String Chromosome;
    public final int Position;
    public final byte Orientation;

    public final String OtherChromosome;
    public final int OtherPosition;
    public final byte OtherOrientation;

    public final int[] Cipos;
    public final int[] Ihompos;
    public final String Homseq;
    public final String InsertSequence;

    public final String SvType;
    public final Set<String> Filters;
    public final VcfType SourceVcfType;

    public VariantBreakend(final VariantContext context, SvCaller svCaller, VcfType sourceVcfType)
    {
        Context = context;

        String alt = context.getAlternateAllele(0).getDisplayString();
        AltCoords = VariantAltInsertCoords.fromRefAlt(alt, alt.substring(0, 1));

        Chromosome = Context.getContig();
        Position = Context.getStart();
        Orientation = AltCoords.Orient.asByte();

        OtherChromosome = AltCoords.OtherChromsome;
        OtherPosition = AltCoords.OtherPosition;
        OtherOrientation = AltCoords.OtherOrient == null ? 0 : AltCoords.OtherOrient.asByte();

        Cipos = getPositionOffsets(CIPOS);
        Ihompos = getPositionOffsets(IHOMPOS);
        Homseq = context.getAttributeAsString(HOMSEQ, "");
        InsertSequence = AltCoords.InsertSequence;

        SvType = svCaller == SvCaller.GRIDSS ?
                context.getAttributeAsString(EVENT_TYPE, "") :
                context.getAttributeAsString(SVTYPE, "");

        Filters = Context.getFilters();

        SourceVcfType = sourceVcfType;
    }

    private boolean isSingle()
    {
        return AltCoords.OtherChromsome.isEmpty();
    }

    private boolean isEnd()
    {
        if(OtherChromosome.isEmpty())
        {
            return false;
        }

        if(Context.getContig().equals(OtherChromosome))
        {
            return Position > OtherPosition;
        }

        return HumanChromosome.lowerChromosome(OtherChromosome, Chromosome);
    }

    public boolean isInverted()
    {
        return Orientation == OtherOrientation;
    }

    private int[] getPositionOffsets(String vcfTag)
    {
        List<Integer> offsetsList = Context.getAttributeAsIntList(vcfTag, 0);
        return offsetsList.size() == 2 ?
                new int[] { offsetsList.get(0), offsetsList.get(1) } :
                new int[] { 0, 0 };
    }

    public int minPosition()
    {
        return Position + Cipos[0];
    }

    public int maxPosition()
    {
        return Position + Cipos[1];
    }

    public int otherMinPosition()
    {
        return isInverted() ?
                OtherPosition - Cipos[0] :
                OtherPosition + Cipos[0];
    }

    public int otherMaxPosition()
    {
        return isInverted() ?
                OtherPosition - Cipos[1] :
                OtherPosition + Cipos[1];
    }

    public boolean exactMatch(final VariantBreakend otherBreakend)
    {
        return
                // First side
                Chromosome.equals(otherBreakend.Chromosome) &&
                        Orientation == otherBreakend.Orientation &&
                        minPosition() == otherBreakend.minPosition() &&
                        maxPosition() == otherBreakend.maxPosition() &&
                        InsertSequence.equals(otherBreakend.InsertSequence) &&
                        Homseq.equals(otherBreakend.Homseq) &&
                        SvType.equals(otherBreakend.SvType) &&

                        // Second side
                        OtherChromosome.equals(otherBreakend.OtherChromosome) &&
                        OtherOrientation == otherBreakend.OtherOrientation &&
                        otherMinPosition() == otherBreakend.otherMinPosition() &&
                        otherMaxPosition() == otherBreakend.otherMaxPosition()
                // No need to check insert seq, hom seq or SV type again. It is always the same on the other side
                ;
    }

    public boolean coordsOnlyMatch(final VariantBreakend otherBreakend)
    {
        return
                // First side
                Chromosome.equals(otherBreakend.Chromosome) &&
                        Position == otherBreakend.Position &&
                        Orientation == otherBreakend.Orientation &&

                        // Second side
                        OtherChromosome.equals(otherBreakend.OtherChromosome) &&
                        OtherPosition == otherBreakend.OtherPosition &&
                        OtherOrientation == otherBreakend.OtherOrientation;
    }

    public boolean approxMatch(final VariantBreakend otherBreakend, final int upperLowerBound)
    {
        return
                // First side
                Chromosome.equals(otherBreakend.Chromosome) &&
                        Orientation == otherBreakend.Orientation &&
                        positionWithin(Position, otherBreakend.Position - upperLowerBound, otherBreakend.Position + upperLowerBound) &&

                        // Second side
                        OtherChromosome.equals(otherBreakend.OtherChromosome) &&
                        OtherOrientation == otherBreakend.OtherOrientation &&
                        positionWithin(OtherPosition,
                                otherBreakend.OtherPosition - upperLowerBound, otherBreakend.OtherPosition + upperLowerBound)
                ;
    }

    public String getExtendedAttributeAsString(String id, String key)
    {
        Object value = Context.getGenotype(id).getExtendedAttribute(key);
        return value != null ? value.toString() : "";
    }

    public int getExtendedAttributeAsInt(String id, String key)
    {
        return getGenotypeAttributeAsInt(Context.getGenotype(id), key, 0);
    }

    public double getExtendedAttributeAsDouble(String id, String key)
    {
        return getGenotypeAttributeAsDouble(Context.getGenotype(id), key, 0);
    }

    public String toString()
    {
        return String.format("coords(%s) cipos(%d,%d)", coordStr(), Cipos[0], Cipos[1]);
    }

    public String coordStr()
    {
        return String.format("%s:%d:%d", Context.getContig(), Position, Orientation);
    }

    public String otherCoordStr()
    {
        return String.format("%s:%d:%d", OtherChromosome, OtherPosition, OtherOrientation);
    }

    public String svCoordStr()
    {
        if(isSingle())
        {
            return coordStr();
        }

        if(isEnd())
        {
            return otherCoordStr() + "_" + coordStr();
        }

        return coordStr() + "_" + otherCoordStr();
    }

    public String filtersStr()
    {
        return isPassVariant() ? PASS : String.join(",", Filters);
    }

    public String qualStr()
    {
        return String.format("%.0f", Context.getPhredScaledQual());
    }

    public StructuralVariantType svType()
    {
        if(OtherChromosome.equals(""))
        {
            return SGL;
        }

        return formSvType(
                Chromosome, OtherChromosome,
                Position, OtherPosition,
                AltCoords.Orient, AltCoords.OtherOrient,
                InsertSequence.isEmpty()
        );
    }

    public boolean isPassVariant()
    {
        return Filters.isEmpty();
    }
}

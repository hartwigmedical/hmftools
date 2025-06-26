package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ALLELE_FRACTION;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LINKS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_REPEAT_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.UNIQUE_FRAG_POSITIONS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsDouble;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.fromRefAlt;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class Breakend
{
    public final String VcfId;
    public final VariantContext Context;
    public final String Chromosome;
    public final int Position;
    public final Orientation Orient;
    public final boolean IsStart; // the start breakend in an SV, or true if a SGL

    public final Genotype RefGenotype;
    public final Genotype TumorGenotype;

    public final String InsertSequence;

    public final Interval ConfidenceInterval;
    public final Interval InexactHomology;
    public final boolean IsLineInsertion;

    private final Variant mVariant;

    private Breakend mLineSiteBreakend;

    public Breakend(
            final Variant variant, final boolean isStart, final VariantContext context, final String chromosome, final int position,
            final Orientation orientation, final Genotype refGenotype, final Genotype tumorGenotype)
    {
        mVariant = variant;
        VcfId = context.getID();
        Context = context;
        Chromosome = chromosome;
        Position = position;
        Orient = orientation;
        IsStart = isStart;

        RefGenotype = refGenotype;
        TumorGenotype = tumorGenotype;

        ConfidenceInterval = Interval.fromCiposTag(context.getAttributeAsIntList(CIPOS, 0));

        String ref = context.getAlleles().get(0).getDisplayString();
        VariantAltInsertCoords altInsertCoords = fromRefAlt(context.getAlleles().get(1).getDisplayString(), ref);
        InsertSequence = altInsertCoords.InsertSequence;

        IsLineInsertion = isMobileLineElement(orientation.asByte(), InsertSequence);

        if(context.hasAttribute(IHOMPOS))
        {
            final List<Integer> ihompos = context.getAttributeAsIntList(IHOMPOS, 0);
            InexactHomology = new Interval(ihompos.get(0), ihompos.get(1));
        }
        else
        {
            InexactHomology = new Interval();
        }

        mLineSiteBreakend = null;
    }

    public static Breakend from(
            final Variant variant, final boolean isStart, final StructuralVariantLeg svLeg,
            final VariantContext variantContext, final int referenceOrdinal, final int tumorOrdinal)
    {
        final Genotype tumorGenotype = tumorOrdinal >= 0 ? variantContext.getGenotype(tumorOrdinal) : null;
        final Genotype refGenotype = referenceOrdinal >= 0 ? variantContext.getGenotype(referenceOrdinal) : null;

        return new Breakend(
                variant, isStart, variantContext, svLeg.chromosome(), svLeg.position(), Orientation.fromByte(svLeg.orientation()),
                refGenotype, tumorGenotype);
    }

    public Variant sv() { return mVariant; }

    public Breakend otherBreakend()
    {
        if(mVariant.isSgl())
            return null;

        return IsStart ? mVariant.breakendEnd() : mVariant.breakendStart();
    }

    public boolean isStart() { return IsStart;}
    public boolean isEnd() { return !IsStart;}

    public double calcAllelicFrequency(final Genotype genotype)
    {
        // set in the depth annotator, which has the same logic as here - so can remove this in future
        if(genotype.hasExtendedAttribute(ALLELE_FRACTION))
            return getGenotypeAttributeAsDouble(genotype, ALLELE_FRACTION, 0);

        int readPairSupport = (mVariant.isSgl() || !mVariant.isShortLocal()) ? getGenotypeAttributeAsInt(genotype, REF_DEPTH_PAIR, 0) : 0;
        int refSupport = getGenotypeAttributeAsInt(genotype, REF_DEPTH, 0);

        int fragmentCount = fragmentCount(genotype);
        double totalSupport = fragmentCount + refSupport + readPairSupport;

        return totalSupport > 0 ? fragmentCount / totalSupport : 0;
    }

    public int fragmentCount(final Genotype genotype)
    {
        return getGenotypeAttributeAsInt(genotype, TOTAL_FRAGS, 0);
    }

    public int fragmentCount() { return fragmentCount(TumorGenotype) + fragmentCount(RefGenotype); }

    // convenience
    public boolean isSgl() { return mVariant.isSgl(); }
    public StructuralVariantType type() { return mVariant.type(); }

    public boolean isLine() { return IsLineInsertion || Context.hasAttribute(LINE_SITE); }
    public void setLineSiteBreakend(final Breakend breakend) { mLineSiteBreakend = breakend; }
    public Breakend lineSiteBreakend() { return mLineSiteBreakend; }

    public boolean inChainedAssembly() { return Context.hasAttribute(ASM_LINKS); }

    public int anchorLength()
    {
        return Context.getAttributeAsInt(SEG_REPEAT_LENGTH, 0);
    }
    public int uniqueFragmentPositions() { return Context.getAttributeAsInt(UNIQUE_FRAG_POSITIONS, 0); }

    public String toString()
    {
        return String.format("%s:%s pos(%s:%d:%d)", VcfId, type(), Chromosome, Position, Orient.asByte());
    }
}

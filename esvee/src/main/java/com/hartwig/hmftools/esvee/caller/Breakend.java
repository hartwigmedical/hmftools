package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASSEMBLY_LINKS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.REF_DEPTH_PAIR;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.fromRefAlt;

import java.util.Collections;
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
    private final List<String> mLinkedAssemblyIds;
    public double mAllelicFrequency;
    private int mChrLocationIndex;

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
        final VariantAltInsertCoords altInsertCoords = fromRefAlt(context.getAlleles().get(1).getDisplayString(), ref);
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

        if(context.hasAttribute(ASSEMBLY_LINKS))
            mLinkedAssemblyIds = context.getAttributeAsStringList(ASSEMBLY_LINKS, "");
        else
            mLinkedAssemblyIds = Collections.emptyList();

        mChrLocationIndex = -1;
    }

    public static Breakend from(
            final Variant variant, final boolean isStart, final StructuralVariantLeg svLeg,
            final VariantContext variantContext, final int referenceOrdinal, final int tumorOrdinal)
    {
        final Genotype tumorGenotype = variantContext.getGenotype(tumorOrdinal);
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

    public boolean isEnd() { return !mVariant.isSgl() && mVariant.breakendEnd() == this;}

    public double calcAllelicFrequency(final Genotype genotype)
    {
        // TODO: just use AF (tag: ALLELE_FRACTION) if present??

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

    // convenience
    public boolean isSgl() { return mVariant.isSgl(); }
    public StructuralVariantType type() { return mVariant.type(); }

    public int minPosition() { return Position + ConfidenceInterval.Start; }
    public int maxPosition() { return Position + ConfidenceInterval.End; }

    public List<String> getAssemblies() { return mLinkedAssemblyIds; }

    public void setChrLocationIndex(int index) { mChrLocationIndex = index; }
    public int chrLocationIndex() { return mChrLocationIndex; }

    public String toString()
    {
        return String.format("%s:%s pos(%s:%d:%d)", VcfId, type(), Chromosome, Position, Orient.asByte());
    }
}

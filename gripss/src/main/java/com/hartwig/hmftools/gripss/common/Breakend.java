package com.hartwig.hmftools.gripss.common;

import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_BAQ;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_BQ;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IHOMPOS;
import static com.hartwig.hmftools.gripss.common.VcfUtils.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.gripss.common.VcfUtils.sglFragmentCount;
import static com.hartwig.hmftools.gripss.common.VariantAltInsertCoords.parseRefAlt;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_CIPOS;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_CIRPOS;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_QUAL;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REF;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_REFPAIR;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_VF;
import static com.hartwig.hmftools.gripss.common.VcfUtils.parseAssemblies;

import java.util.List;

import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class Breakend
{
    public final String VcfId;
    public final VariantContext Context;
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;
    public final boolean IsStart; // the start breakend in an SV, or true if a SGL

    public final double Qual;
    public final int TumorFragments;
    public final int ReferenceFragments;

    public final Genotype RefGenotype;
    public final Genotype TumorGenotype;

    public final String Ref;
    public final String Alt;
    public final String InsertSequence;
    public final String OtherChromosome;
    public final int OtherPosition;
    public final byte OtherOrientation;

    public final Interval ConfidenceInterval;
    public final Interval RemoteConfidenceInterval;
    public final Interval InexactHomology;
    public final boolean IsLineInsertion;

    private final SvData mSvData;
    private final List<String> mAssemblies;
    private boolean mReligned;
    public double mAllelicFrequency;
    private int mChrLocationIndex;

    public Breakend(
            final SvData svData, final boolean isStart, final VariantContext context, final String chromosome, final int position,
            final byte orientation, final Genotype refGenotype, final Genotype tumorGenotype)
    {
        mSvData = svData;
        VcfId = context.getID();
        Context = context;
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        IsStart = isStart;

        RefGenotype = refGenotype;
        TumorGenotype = tumorGenotype;

        // NOTE: SvData is only partially constructed so be careful which fields are used
        boolean isSgl = mSvData.type() == SGL;

        if(refGenotype != null)
        {
            ReferenceFragments = isSgl ? sglFragmentCount(refGenotype) : getGenotypeAttributeAsInt(refGenotype, VT_VF, 0);
        }
        else
        {
            ReferenceFragments = 0;
        }

        TumorFragments = isSgl ? sglFragmentCount(tumorGenotype) : getGenotypeAttributeAsInt(tumorGenotype, VT_VF, 0);

        setAllelicFrequency();

        ConfidenceInterval = VcfUtils.confidenceInterval(context, VT_CIPOS);
        RemoteConfidenceInterval = VcfUtils.confidenceInterval(context, VT_CIRPOS);

        Ref = context.getAlleles().get(0).getDisplayString();

        final VariantAltInsertCoords altInsertCoords = parseRefAlt(context.getAlleles().get(1).getDisplayString(), Ref);
        Alt = altInsertCoords.Alt;

        InsertSequence = altInsertCoords.InsertSequence;
        OtherChromosome = altInsertCoords.Chromsome;
        OtherPosition = altInsertCoords.Position;
        OtherOrientation = altInsertCoords.Orientation;

        IsLineInsertion = isMobileLineElement(orientation, InsertSequence);

        if(mSvData.type() == SGL)
        {
            final String qualTag = IsLineInsertion ? VT_BQ : VT_BAQ;
            Qual = VcfUtils.getGenotypeAttributeAsDouble(tumorGenotype, qualTag, 0);
        }
        else
        {
            Qual = VcfUtils.getGenotypeAttributeAsDouble(tumorGenotype, VT_QUAL, 0);
        }

        if(context.hasAttribute(VT_IHOMPOS))
        {
            final List<Integer> ihompos = context.getAttributeAsIntList(VT_IHOMPOS, 0);
            InexactHomology = new Interval(ihompos.get(0), ihompos.get(1));
        }
        else
        {
            InexactHomology = new Interval();
        }

        mAssemblies = parseAssemblies(context);
        mReligned = false;
        mChrLocationIndex = -1;
    }

    public static Breakend from(
            final SvData svData, final boolean isStart, final StructuralVariantLeg svLeg,
            final VariantContext variantContext, final int referenceOrdinal, final int tumorOrdinal)
    {
        final Genotype tumorGenotype = variantContext.getGenotype(tumorOrdinal);
        final Genotype refGenotype = referenceOrdinal >= 0 ? variantContext.getGenotype(referenceOrdinal) : null;

        return new Breakend(
                svData, isStart, variantContext, svLeg.chromosome(), (int)svLeg.position(), svLeg.orientation(), refGenotype, tumorGenotype);
    }

    public static Breakend realigned(final Breakend original, final VariantContext newContext, final int newPosition)
    {
        Breakend newBreakend = new Breakend(
                original.sv(), original.IsStart, newContext, original.Chromosome, newPosition, original.Orientation,
                original.RefGenotype, original.TumorGenotype);

        newBreakend.markRealigned();
        return newBreakend;
    }

    public final SvData sv() { return mSvData; }

    public Breakend otherBreakend()
    {
        if(mSvData.isSgl())
            return null;

        return IsStart ? mSvData.breakendEnd() : mSvData.breakendStart();
    }

    public double allelicFrequency() { return mAllelicFrequency; }

    public void setAllelicFrequency()
    {
        int readPairSupport = (mSvData.isSgl() || !mSvData.isShortLocal()) ? getGenotypeAttributeAsInt(TumorGenotype, VT_REFPAIR, 0) : 0;
        int refSupport = getGenotypeAttributeAsInt(TumorGenotype, VT_REF, 0);
        double totalSupport = TumorFragments + refSupport + readPairSupport;
        mAllelicFrequency = totalSupport > 0 ? TumorFragments / totalSupport : 0;
    }

    // convenience
    public boolean isSgl() { return mSvData.isSgl(); }
    public StructuralVariantType type() { return mSvData.type(); }
    public boolean imprecise() { return mSvData.imprecise(); }
    public boolean posOrient() { return Orientation == POS_ORIENT; }

    public int insertSequenceLength() { return InsertSequence.length(); }

    public int minPosition() { return Position + ConfidenceInterval.Start; }
    public int maxPosition() { return Position + ConfidenceInterval.End; }

    public List<String> getAssemblies() { return mAssemblies; }

    public void markRealigned() { mReligned = true; }
    public boolean realigned() { return mReligned; }

    public void setChrLocationIndex(int index) { mChrLocationIndex = index; }
    public int chrLocationIndex() { return mChrLocationIndex; }

    public String toString()
    {
        return String.format("%s:%s pos(%s:%d:%d)", VcfId, type(), Chromosome, Position, Orientation);
    }
}

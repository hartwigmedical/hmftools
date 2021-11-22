package com.hartwig.hmftools.gripss.common;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BAQ;
import static com.hartwig.hmftools.gripss.VcfUtils.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.gripss.VcfUtils.sglFragmentCount;
import static com.hartwig.hmftools.gripss.common.VariantAltInsertCoords.parseRefAlt;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BQ;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_BVF;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_CIPOS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_CIRPOS;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_QUAL;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_REF;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_REFPAIR;
import static com.hartwig.hmftools.gripss.VcfUtils.VT_VF;
import static com.hartwig.hmftools.gripss.VcfUtils.parseAssemblies;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.gripss.VcfUtils;
import com.hartwig.hmftools.gripss.filters.FilterType;

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
    public final int ReferenceReads;
    public final int ReferencePairReads;

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

    private final SvData mSvData;
    private final List<FilterType> mFilters;
    private final List<String> mAssemblies;
    private boolean mReligned;
    private int mChrLocationIndex;

    public Breakend(
            final SvData svData, final boolean isStart, final VariantContext context, final String chromosome, final int position, final byte orientation,
            final Genotype refGenotype, final Genotype tumorGenotype, final double qual, final int tumorFragments,
            final int refFrags, final int refReads, final int refPairReads)
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
        Qual = qual;
        TumorFragments = tumorFragments;
        ReferenceFragments = refFrags;
        ReferenceReads = refReads;
        ReferencePairReads = refPairReads;

        ConfidenceInterval = VcfUtils.confidenceInterval(context, VT_CIPOS);
        RemoteConfidenceInterval = VcfUtils.confidenceInterval(context, VT_CIRPOS);

        Ref = context.getAlleles().get(0).getDisplayString();;

        final VariantAltInsertCoords altInsertCoords = parseRefAlt(context.getAlleles().get(1).getDisplayString(), Ref);
        Alt = altInsertCoords.Alt;

        InsertSequence = altInsertCoords.InsertSequence;
        OtherChromosome = altInsertCoords.Chromsome;
        OtherPosition = altInsertCoords.Position;
        OtherOrientation = altInsertCoords.Orientation;

        mFilters = Lists.newArrayList();
        mAssemblies = parseAssemblies(context);
        mReligned = false;
        mChrLocationIndex = -1;
    }

    public static Breakend from(
            final SvData svData, final StructuralVariantType type, final boolean isStart, final StructuralVariantLeg svLeg,
            final VariantContext variantContext, final int referenceOrdinal, final int tumorOrdinal)
    {
        final Genotype tumorGenotype = variantContext.getGenotype(tumorOrdinal);
        final Genotype refGenotype = referenceOrdinal >= 0 ? variantContext.getGenotype(referenceOrdinal) : null;

        final String qualTag = type == SGL ? VT_BAQ : VT_QUAL;
        final String fragsTag = type == SGL ? VT_BVF : VT_VF;
        double qual = VcfUtils.getGenotypeAttributeAsDouble(tumorGenotype, qualTag, 0);

        int refFrags = 0;
        int refReads = 0;
        int refPairReads = 0;

        if(refGenotype != null)
        {
            refFrags = type == SGL ? sglFragmentCount(refGenotype) : getGenotypeAttributeAsInt(refGenotype, fragsTag, 0);
            refReads = getGenotypeAttributeAsInt(refGenotype, VT_REF, 0);
            refPairReads = getGenotypeAttributeAsInt(refGenotype, VT_REFPAIR, 0);
        }

        int tumorFrags = type == SGL ? sglFragmentCount(tumorGenotype) : getGenotypeAttributeAsInt(tumorGenotype, fragsTag, 0);

        return new Breakend(
                svData, isStart, variantContext, svLeg.chromosome(), (int)svLeg.position(), svLeg.orientation(), refGenotype, tumorGenotype,
                qual, tumorFrags, refFrags, refReads, refPairReads);
    }

    public static Breakend realigned(final Breakend original, final VariantContext newContext, final int newPosition)
    {
        Breakend newBreakend = new Breakend(
                original.sv(), original.IsStart, newContext, original.Chromosome, newPosition, original.Orientation, original.RefGenotype,
                original.TumorGenotype, original.Qual, original.TumorFragments, original.ReferenceFragments, original.ReferenceReads,
                original.ReferencePairReads);

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

    // convenience
    public boolean isSgl() { return mSvData.isSgl(); }
    public StructuralVariantType type() { return mSvData.type(); }
    public String svId() { return mSvData.id(); }
    public boolean imprecise() { return mSvData.imprecise(); }
    public boolean posOrient() { return Orientation == POS_ORIENT; }
    public boolean negOrient() { return Orientation == NEG_ORIENT; }

    public int insertSequenceLength() { return InsertSequence.length(); }

    public int minPosition() { return Position + ConfidenceInterval.Start; }
    public int maxPosition() { return Position + ConfidenceInterval.End; }

    public List<FilterType> getFilters() { return mFilters; }

    public void addFilter(final FilterType filter)
    {
        if(!mFilters.contains(filter))
            mFilters.add(filter);
    }

    public List<String> getAssemblies() { return mAssemblies; }

    public void markRealigned() { mReligned = true; }
    public boolean realigned() { return mReligned; }

    public void setChrLocationIndex(int index) { mChrLocationIndex = index; }
    public int chrLocationIndex() { return mChrLocationIndex; }

    /*
    val startBreakend: Breakend = Breakend(contig, start + confidenceInterval.first, start + confidenceInterval.second, orientation)
    val endBreakend: Breakend? = (variantType as? Paired)?.let { Breakend(it.otherChromosome, it.otherPosition + remoteConfidenceInterval.first, it.otherPosition + remoteConfidenceInterval.second, it.endOrientation) }
     */


}

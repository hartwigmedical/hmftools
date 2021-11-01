package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.CIPOS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.svtools.germline.VariantAltInsertCoords.parseRefAlt;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BQ;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.BVF;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.CIRPOS;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.QUAL;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.REF;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.REFPAIR;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.VF;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.parseAssemblies;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class Breakend
{
    public final String SvId;
    public final String VcfId;
    public final VariantContext Context;
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;

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

    private final List<FilterType> mFilters;
    private final List<String> mAssemblies;
    private boolean mReligned;

    public Breakend(
            final String svId, final VariantContext context, final String chromosome, final int position, final byte orientation,
            final Genotype refGenotype, final Genotype tumorGenotype, final double qual, final int tumorFragments,
            final int refFrags, final int refReads, final int refPairReads)
    {
        VcfId = context.getID();
        SvId = svId;
        Context = context;
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;

        RefGenotype = refGenotype;
        TumorGenotype = tumorGenotype;
        Qual = qual;
        TumorFragments = tumorFragments;
        ReferenceFragments = refFrags;
        ReferenceReads = refReads;
        ReferencePairReads = refPairReads;

        ConfidenceInterval = VcfUtils.confidenceInterval(context, CIPOS);
        RemoteConfidenceInterval = VcfUtils.confidenceInterval(context, CIRPOS);

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
    }

    public static Breakend from(
            final String svId, final StructuralVariantType type, final StructuralVariantLeg svLeg, final VariantContext variantContext,
            final int referenceOrdinal, final int tumorOrdinal)
    {
        final Genotype tumorGenotype = variantContext.getGenotype(tumorOrdinal);
        final Genotype refGenotype = referenceOrdinal > 0 ? variantContext.getGenotype(referenceOrdinal) : null;

        final String qualTag = type == SGL ? BQ : QUAL;
        final String fragsTag = type == SGL ? BVF : VF;
        double qual = VcfUtils.getGenotypeAttributeAsDouble(tumorGenotype, qualTag, 0);

        int refFrags = 0;
        int refReads = 0;
        int refPairReads = 0;

        if(refGenotype != null)
        {
            refFrags = VcfUtils.getGenotypeAttributeAsInt(refGenotype, fragsTag, 0);
            refReads = VcfUtils.getGenotypeAttributeAsInt(refGenotype, REF, 0);
            refPairReads = VcfUtils.getGenotypeAttributeAsInt(refGenotype, REFPAIR, 0);
        }

        int tumorFrags = VcfUtils.getGenotypeAttributeAsInt(tumorGenotype, fragsTag, 0);

        return new Breakend(
                svId, variantContext, svLeg.chromosome(), (int)svLeg.position(), svLeg.orientation(), refGenotype, tumorGenotype,
                qual, tumorFrags, refFrags, refReads, refPairReads);
    }

    public static Breakend realigned(final Breakend original, final VariantContext newContext, final int newPosition)
    {
        Breakend newBreakend = new Breakend(
                original.VcfId, newContext, original.Chromosome, newPosition, original.Orientation, original.RefGenotype,
                original.TumorGenotype, original.Qual, original.TumorFragments, original.ReferenceFragments, original.ReferenceReads,
                original.ReferencePairReads);

        newBreakend.markRealigned();
        return newBreakend;
    }

    public int insertSequenceLength() { return 0; } // TODO

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

    /*
    val startBreakend: Breakend = Breakend(contig, start + confidenceInterval.first, start + confidenceInterval.second, orientation)
    val endBreakend: Breakend? = (variantType as? Paired)?.let { Breakend(it.otherChromosome, it.otherPosition + remoteConfidenceInterval.first, it.otherPosition + remoteConfidenceInterval.second, it.endOrientation) }
     */


}

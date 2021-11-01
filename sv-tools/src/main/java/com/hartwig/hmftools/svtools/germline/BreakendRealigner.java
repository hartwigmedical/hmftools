package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.CIPOS;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.svtools.germline.VariantAltInsertCoords.formPairedAltString;
import static com.hartwig.hmftools.svtools.germline.VariantAltInsertCoords.formSingleAltString;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.CIRPOS;
import static com.hartwig.hmftools.svtools.germline.VcfUtils.REALIGN;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class BreakendRealigner
{
    private final RefGenomeInterface mRefGenome;
    private final int mReferenceOrdinal;
    private final int mTumorOrdinal;

    // comparator: ContigComparator

    public BreakendRealigner(final RefGenomeInterface refGenome, final int referenceOrdinal, final int tumorOrdinal)
    {
        mRefGenome = refGenome;
        mReferenceOrdinal = referenceOrdinal;
        mTumorOrdinal = tumorOrdinal;
    }

    public Breakend realignRemote(final Breakend breakend, final Breakend realignedOther)
    {
        Allele refAllele = breakend.Context.getReference();

        List<Allele> alleles = Lists.newArrayList();
        alleles.add(refAllele);

        // TODO
        // val mate = variantType as Paired
        // String newAlt = ""; // mate.altString(other.Position, refAllele.displayString)

        String newAlt = formPairedAltString(
                refAllele.getDisplayString(), breakend.InsertSequence,
                breakend.OtherChromosome, realignedOther.Position, breakend.Orientation, realignedOther.Orientation);

        alleles.add(Allele.create(newAlt));

        VariantContext newContext = new VariantContextBuilder(breakend.Context)
                .alleles(alleles)
                .attribute(CIRPOS, Lists.newArrayList(breakend.ConfidenceInterval.Start, breakend.ConfidenceInterval.End))
                .make();

        return Breakend.realigned(breakend, newContext, breakend.Position);
    }

    public Breakend realign(final Breakend breakend, boolean isSgl, boolean imprecise)
    {
        if(breakend.insertSequenceLength() > 0)
            return breakend;

        if(!imprecise && !isSgl)
        {
            Interval[] intervals = centreAlignConfidenceIntervals(breakend);
            Interval centeredCipos = intervals[0];
            Interval centeredRemoteCipos = intervals[1];

            if(!centeredCipos.matches(breakend.ConfidenceInterval) || !centeredRemoteCipos.matches(breakend.RemoteConfidenceInterval))
            {
                return realignPaired(breakend, centeredCipos);
            }
        }

        if(imprecise && (breakend.ConfidenceInterval.Start != 0 || breakend.ConfidenceInterval.End != 0))
        {
            return isSgl ? sideAlignSingle(breakend) : sideAlignPaired(breakend);
        }

        return breakend;
    }

    private Interval[] centreAlignConfidenceIntervals(final Breakend breakend)
    {
        // val mate = variantType as Paired
        final Breakend otherBreakend = null;

        // TODO
        boolean invertStart = breakend.Orientation == otherBreakend.Orientation;
        // && comparator.compare(contig, start, mate.otherChromosome, mate.otherPosition) > 0

        boolean invertEnd = breakend.Orientation == otherBreakend.Orientation && !invertStart;

        Interval centeredCipos = centreAlignConfidenceInterval(invertStart, breakend.ConfidenceInterval);
        Interval centeredRemoteCipos = centreAlignConfidenceInterval(invertEnd, breakend.RemoteConfidenceInterval);
        return new Interval[] { centeredCipos, centeredRemoteCipos };
    }

    private Breakend sideAlignPaired(final Breakend breakend)
    {
        Interval newCipos = sideAlignConfidenceInterval(breakend.Orientation, breakend.ConfidenceInterval);
        return realignPaired(breakend, newCipos);
    }

    private Breakend realignPaired(final Breakend breakend, final Interval newCipos)
    {
        if(newCipos.matches(breakend.ConfidenceInterval))
            return breakend;

        int newStart = updatedPosition(breakend.Position, breakend.ConfidenceInterval, newCipos);

        String newRef = mRefGenome.getBaseString(breakend.Chromosome, newStart, newStart);
        // getSubsequenceAt(contig, newStart.toLong(), newStart.toLong()).unambiguousNucleotides

        // val mate = variantType as Paired
        // String newAlt = ""; // mate.altString(newRef);

        String newAlt = formPairedAltString(
                newRef, breakend.InsertSequence, breakend.OtherChromosome, breakend.Position, breakend.Orientation, breakend.OtherOrientation);

        List<Allele> alleles = Lists.newArrayList();
        alleles.add(Allele.create(newRef, true));
        alleles.add(Allele.create(newAlt));

        VariantContext newContext = new VariantContextBuilder(breakend.Context)
                .start(newStart)
                .stop(newStart)
                .alleles(alleles)
                .attribute(REALIGN, true)
                .attribute(CIPOS, Lists.newArrayList(newCipos.Start, newCipos.End))
                .make();

        return Breakend.realigned(breakend, newContext, newStart);
    }

    private Breakend sideAlignSingle(final Breakend breakend)
    {
        Interval newCipos = sideAlignConfidenceInterval(breakend.Orientation, breakend.ConfidenceInterval);
        int newStart = updatedPosition(breakend.Position, breakend.ConfidenceInterval, newCipos);

        String newRef = mRefGenome.getBaseString(breakend.Chromosome, newStart, newStart);
        // getSubsequenceAt(contig, newStart.toLong(), newStart.toLong()).unambiguousNucleotides
        // TODO: understand unambiguousNucleotides

        // val mate = variantType as Single
        // String newAlt = ""; // mate.altString(newRef);
        String newAlt = formSingleAltString(newRef, breakend.InsertSequence, breakend.Orientation);
        // val alleles = listOf(Allele.create(newRef, true), Allele.create(mate.altString(newRef)))

        List<Allele> alleles = Lists.newArrayList();
        alleles.add(Allele.create(newRef, true));
        alleles.add(Allele.create(newAlt));

        VariantContext newContext = new VariantContextBuilder(breakend.Context)
                .start(newStart)
                .stop(newStart)
                .alleles(alleles)
                .attribute(REALIGN, true)
                .attribute(CIPOS, Lists.newArrayList(newCipos.Start, newCipos.End))
                .make();

        return Breakend.realigned(breakend, newContext, newStart);
    }

    private static int updatedPosition(int position, final Interval oldCipos, final Interval newCipos)
    {
        return position + oldCipos.Start - newCipos.Start;
    }

    private static Interval centreAlignConfidenceInterval(boolean invert, final Interval cipos)
    {
        int totalRange = cipos.End - cipos.Start;
        int newCiposStart = -totalRange / 2;
        int newCiposEnd = totalRange + newCiposStart;

        if(invert)
            return new Interval(-newCiposEnd, -newCiposStart);
        else
            return new Interval(newCiposStart, newCiposEnd);
    }

    private static Interval sideAlignConfidenceInterval(byte orientation, final Interval cipos)
    {
        if(orientation == POS_ORIENT)
            return new Interval(0, cipos.End - cipos.Start);
        else
            return new Interval(cipos.Start - cipos.End, 0);
    }
}

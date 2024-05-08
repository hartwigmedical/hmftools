package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome.chromosomeRank;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formPairedAltString;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.formSingleAltString;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.CIRPOS;
import static com.hartwig.hmftools.common.sv.gridss.GridssVcfTags.REALIGN;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.Interval;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class BreakendRealigner
{
    private final RefGenomeInterface mRefGenome;

    public BreakendRealigner(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
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

    public Breakend realignRemote(final Breakend breakend, final Breakend realignedOther)
    {
        // realign the breakend, using the already realigned other breakend's new coords
        Allele refAllele = breakend.Context.getReference();

        List<Allele> alleles = Lists.newArrayList();
        alleles.add(refAllele);

        String newAlt = formPairedAltString(
                refAllele.getDisplayString(), breakend.InsertSequence,
                realignedOther.Chromosome, realignedOther.Position,
                Orientation.fromByte(breakend.Orientation), Orientation.fromByte(realignedOther.Orientation));

        alleles.add(Allele.create(newAlt));

        VariantContext newContext = new VariantContextBuilder(breakend.Context)
                .alleles(alleles)
                .attribute(CIRPOS, Lists.newArrayList(realignedOther.ConfidenceInterval.Start, realignedOther.ConfidenceInterval.End))
                .make();

        return Breakend.realigned(breakend, newContext, breakend.Position);
    }

    private Interval[] centreAlignConfidenceIntervals(final Breakend breakend)
    {
        // val mate = variantType as Paired
        final Breakend otherBreakend = breakend.otherBreakend();

        boolean invertStart;

        if(breakend.Orientation == otherBreakend.Orientation)
        {
            // comparator: -ve if chromosome is lower or if positon is lower
            if(breakend.Chromosome.equals(otherBreakend.Chromosome))
                invertStart = breakend.Position > otherBreakend.Position;
            else
                invertStart = chromosomeRank(breakend.Chromosome) > chromosomeRank(otherBreakend.Chromosome);
        }
        else
        {
            invertStart = false;
        }

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

        String newRef = toStandardNucleotides(mRefGenome.getBaseString(breakend.Chromosome, newStart, newStart));

        final Breakend otherBreakend = breakend.otherBreakend();

        String newAlt = formPairedAltString(
                newRef, otherBreakend.InsertSequence, otherBreakend.Chromosome, otherBreakend.Position,
                Orientation.fromByte(breakend.Orientation), Orientation.fromByte(otherBreakend.Orientation));

        List<Allele> alleles = Lists.newArrayList();
        alleles.add(Allele.create(newRef, true));
        alleles.add(Allele.create(newAlt));

        VariantContextBuilder builder = new VariantContextBuilder(breakend.Context)
                .start(newStart)
                .stop(newStart)
                .alleles(alleles)
                .attribute(REALIGN, true)
                .attribute(CIPOS, Lists.newArrayList(newCipos.Start, newCipos.End));

        if(breakend.InexactHomology.length() > 0)
        {
            int ciposShift = newCipos.Start - breakend.ConfidenceInterval.Start;

            builder.attribute(IHOMPOS, Lists.newArrayList(
                    breakend.InexactHomology.Start + ciposShift, breakend.InexactHomology.End + ciposShift));
        }

        VariantContext newContext = builder.make();

        return Breakend.realigned(breakend, newContext, newStart);
    }

    private Breakend sideAlignSingle(final Breakend breakend)
    {
        Interval newCipos = sideAlignConfidenceInterval(breakend.Orientation, breakend.ConfidenceInterval);
        int newStart = updatedPosition(breakend.Position, breakend.ConfidenceInterval, newCipos);

        String newRef = toStandardNucleotides(mRefGenome.getBaseString(breakend.Chromosome, newStart, newStart));

        String newAlt = formSingleAltString(newRef, breakend.InsertSequence, Orientation.fromByte(breakend.Orientation));

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

    private static final List<Character> VALID_BASES = Lists.newArrayList('A', 'G', 'C', 'T', 'N');

    private static String toStandardNucleotides(final String bases)
    {
        StringBuilder newBases = new StringBuilder(bases.length());

        for(int i = 0; i < bases.length(); ++i)
        {
            char base = bases.charAt(i);
            if(VALID_BASES.contains(base))
                newBases.append(base);
            else
                newBases.append('N');
        }

        return newBases.toString();
    }
}

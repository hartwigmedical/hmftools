package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.common.CommonUtils.compareJunctions;
import static com.hartwig.hmftools.esvee.common.CommonUtils.formSvType;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.filters.FilterType;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;

public class Breakend implements Comparable<Breakend>
{
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;
    public final String InsertedBases;
    public final HomologyData Homology;

    private Breakend mOtherBreakend;

    private List<AlignData> mAlternativeAlignments;
    private List<BreakendSegment> mSegments;

    private final List<BreakendSupport> mBreakendSupport; // one for each sample loaded, indexed as per config
    private final Set<FilterType> mFilters;

    public Breakend(
            final String chromosome, final int position, final byte orientation, final String insertedBases, final HomologyData homology)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        InsertedBases = insertedBases;
        Homology = homology;

        mOtherBreakend = null;

        mBreakendSupport = Lists.newArrayList();
        mSegments = Lists.newArrayList();

        mAlternativeAlignments = null;

        mFilters = Sets.newHashSet();
    }

    public Breakend otherBreakend() { return mOtherBreakend; }
    public void setOtherBreakend(final Breakend breakend) { mOtherBreakend = breakend; }

    public List<BreakendSegment> segments() { return mSegments; }
    public void addSegment(final BreakendSegment segment) { mSegments.add(segment); }

    public List<AlignData> alternativeAlignments() { return mAlternativeAlignments; }
    public void setAlternativeAlignments(final List<AlignData> altAlignments) { mAlternativeAlignments = altAlignments; }

    public List<BreakendSupport> sampleSupport() { return mBreakendSupport; }

    public StructuralVariantType svType()
    {
        if(mOtherBreakend == null)
            return SGL;

        return formSvType(
                Chromosome, mOtherBreakend.Chromosome, Position, mOtherBreakend.Position,
                Orientation, mOtherBreakend.Orientation, !InsertedBases.isEmpty());
    }

    public int calcSvQual()
    {
        int breakendQual = calcQual();
        return mOtherBreakend != null ? breakendQual + mOtherBreakend.calcQual() : breakendQual;
    }

    public int calcQual()
    {
        return mSegments.stream().mapToInt(x -> x.calcQual()).max().orElse(0);
    }

    public Set<FilterType> filters() { return mFilters; }
    public boolean passing() { return mFilters.isEmpty(); }
    public void addFilter(final FilterType filterType) { mFilters.add(filterType); }

    public boolean readSpansJunction(final SupportRead read, boolean isAlternativeAlignment)
    {
        if(!isAlternativeAlignment)
        {
            // first an aligned junction read
            if(read.unclippedStart() < Position && read.unclippedEnd() > Position)
                return true;

            return false;

            /*
            // check if discordant
            if(Orientation == POS_ORIENT && read.alignmentEnd() <= Position)
                return true;
            else if(Orientation == NEG_ORIENT && read.alignmentStart() >= Position)
                return true;
            else
                return false;
            */
        }
        else
        {
            // next a misaligned junction read - crossing the segment boundary
            int readSeqStartIndex = read.linkedAssemblyIndex();
            int readSeqEndIndex = read.linkedAssemblyIndex() + read.baseLength() - 1;

            // look for a read crossing a segment boundary
            for(BreakendSegment segment : mSegments)
            {
                int segmentIndexStart = segment.Alignment.adjustedSequenceStart();
                int segmentIndexEnd = segment.Alignment.adjustedSequenceEnd();

                if(segmentIndexStart > 0 && readSeqStartIndex < segmentIndexStart && readSeqEndIndex > segmentIndexStart)
                    return true;

                if(segmentIndexEnd < segment.SequenceLength - 1 && readSeqStartIndex < segmentIndexEnd && readSeqEndIndex > segmentIndexEnd)
                    return true;
            }

            return false;
        }
    }

    public boolean matches(final String chromosome, final int position, final byte orientation)
    {
        if(!Chromosome.equals(chromosome) || Orientation != orientation)
            return false;

        if(Homology != null)
            return positionWithin(position, Position + Homology.InexactStart, Position + Homology.InexactEnd);
        else
            return Position == position;
    }

    public String toString()
    {
        return format("%s:%d:%d %s", Chromosome, Position, Orientation, svType());
    }

    @Override
    public int compareTo(final Breakend other)
    {
        return compareJunctions(Chromosome, other.Chromosome, Position, other.Position, Orientation, other.Orientation);
    }
}

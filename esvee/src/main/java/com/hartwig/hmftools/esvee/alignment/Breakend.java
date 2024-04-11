package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.esvee.common.CommonUtils.compareJunctions;
import static com.hartwig.hmftools.esvee.common.CommonUtils.formSvType;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.esvee.assembly.filters.FilterType;

public class Breakend implements Comparable<Breakend>
{
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;
    public final String InsertedBases;
    public final HomologyData Homology;

    private Breakend mOtherBreakend;

    private int mAnchorLength;
    private List<AlignData> mAlternativeAlignments;

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

        mAlternativeAlignments = null;
        mAnchorLength = 0;

        mFilters = Sets.newHashSet();
    }

    public Breakend otherBreakend() { return mOtherBreakend; }
    public void setOtherBreakend(final Breakend breakend) { mOtherBreakend = breakend; }

    public int anchorLength() { return mAnchorLength; }
    public void setAnchorLength(int length) { mAnchorLength = length; }

    public List<AlignData> alternativeAlignments() { return mAlternativeAlignments; }
    public void setAlternativeAlignments(final List<AlignData> altAlignments) { mAlternativeAlignments = altAlignments; }

    public void addBreakendSupport(final BreakendSupport breakendSupport) { mBreakendSupport.add(breakendSupport); }
    public List<BreakendSupport> sampleSupport() { return mBreakendSupport; }

    public StructuralVariantType svType()
    {
        if(mOtherBreakend == null)
            return SGL;

        return formSvType(
                Chromosome, mOtherBreakend.Chromosome, Position, mOtherBreakend.Position,
                Orientation, mOtherBreakend.Orientation, !InsertedBases.isEmpty());
    }

    public Set<FilterType> filters() { return mFilters; }
    public boolean passing() { return mFilters.isEmpty(); }
    public void addFilter(final FilterType filterType) { mFilters.add(filterType); }

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

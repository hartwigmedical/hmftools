package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import java.util.List;

import com.hartwig.hmftools.common.sv.StructuralVariantType;

public class Breakend
{
    public final String Chromosome;
    public final int Position;
    public final byte Orientation;
    public final String InsertedBases;
    public final StructuralVariantType Type;

    private Breakend mLinkedBreakend;

    // other annotations for VCF
    private final String mHomology;
    private final int[] mConfidenceInterval;
    private final int[] mInexactHomology;

    //

    /*
    private final List<Integer> mAssemblyIds;
    private int mAssemblyLength;
    private int mAlignmentLength;
    private int
    Unique id(s) of assembly(s) containing the breakend
    Total length(s) of assembly(s) containing the breakend
    #(s) of segments in assembly(s) containing the breakend
    Breakend position(s) in assembly(s)
    Breakend orientation(s) in assembly(s)
    Unique id(s) of segment(s) containing the breakend
    Aligned length of segment in reference genome
    Highest MAPQ of local aligments
    Highest alignment score of local alignments
    Repeat length of local alignment with highest alignment score
     */

    public Breakend(
            final String chromosome, final int position, final byte orientation, final String insertedBases,
            final StructuralVariantType type)
    {
        Chromosome = chromosome;
        Position = position;
        Orientation = orientation;
        InsertedBases = insertedBases;
        Type = type;

        mHomology = "";
        mConfidenceInterval = new int[2];
        mInexactHomology = new int[2];
        mLinkedBreakend = null;
    }

    public void setOtherBreakend(final Breakend breakend) { mLinkedBreakend = breakend; }


    public String toString()
    {
        return format("%s:%d%d %s", Chromosome, Position, Orientation, Type);
    }
}

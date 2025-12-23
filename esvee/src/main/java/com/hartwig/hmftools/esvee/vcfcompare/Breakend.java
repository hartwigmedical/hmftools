package com.hartwig.hmftools.esvee.vcfcompare;

import static com.hartwig.hmftools.common.sv.LineElements.isMobileLineElement;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_LINKS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.CIPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.IHOMPOS;
import static com.hartwig.hmftools.common.sv.SvVcfTags.LINE_SITE;
import static com.hartwig.hmftools.common.sv.SvVcfTags.PON_FILTER_PON;
import static com.hartwig.hmftools.common.sv.SvVcfTags.SEG_REPEAT_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.VariantAltInsertCoords.fromRefAlt;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.getGenotypeAttributeAsInt;
import static com.hartwig.hmftools.esvee.vcfcompare.CoordMatchType.NONE;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.sv.VariantAltInsertCoords;
import com.hartwig.hmftools.esvee.caller.Interval;

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

    public final Set<String> Filters;

    private final Variant mVariant;

    private CoordMatchType mCoordMatchType;

    private LineData mLineData;

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

        IsLineInsertion = Context.hasAttribute(LINE_SITE);

        if(context.hasAttribute(IHOMPOS))
        {
            List<Integer> ihompos = context.getAttributeAsIntList(IHOMPOS, 0);
            InexactHomology = new Interval(ihompos.get(0), ihompos.size() == 2 ? ihompos.get(1) : 0);
        }
        else
        {
            InexactHomology = new Interval();
        }

        Filters = Sets.newHashSet();

        for(String filterStr : context.getFilters())
        {
            if(!filterStr.equals(PASS))
                Filters.add(filterStr);
        }

        mCoordMatchType = NONE;
        mLineData = null;
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

    public int minPosition()
    {
        return Position + ConfidenceInterval.Start;
    }
    public int maxPosition()
    {
        return Position + ConfidenceInterval.End;
    }

    public int fragmentCount(final Genotype genotype)
    {
        return getGenotypeAttributeAsInt(genotype, TOTAL_FRAGS, 0);
    }

    public int tumorFragmentCount() { return fragmentCount(TumorGenotype); }
    public int referenceFragmentCount() { return fragmentCount(RefGenotype); }

    // convenience
    public boolean isSgl() { return mVariant.isSgl(); }
    public StructuralVariantType type() { return mVariant.type(); }

    public boolean isLine() { return Context.hasAttribute(LINE_SITE); }

    public boolean inChainedAssembly() { return Context.hasAttribute(ASM_LINKS); }

    public int anchorLength()
    {
        return Context.getAttributeAsInt(SEG_REPEAT_LENGTH, 0);
    }

    public boolean isPass() { return Filters.isEmpty(); }
    public boolean isFiltered() { return !isPass() && !isPonOnly(); }
    public boolean isPonOnly() { return Filters.size() == 1 && Filters.contains(PON_FILTER_PON); }

    public void setCoordMatchType(final CoordMatchType type) { mCoordMatchType = type; }
    public CoordMatchType coordMatchType() { return mCoordMatchType; }
    public boolean matched() { return mCoordMatchType != NONE; }

    public void setLineData(final LineData data) { mLineData = data; }
    public LineData lineData() { return mLineData; }
    public boolean hasPolyATail() { return mLineData != null && mLineData.HasPolyAT; }
    public boolean hasLineLink() { return mLineData != null && mLineData.hasLineLink(); }
    public boolean hasInferredLineLink() { return mLineData != null && mLineData.hasInferredLineLink(); }
    public void setLineLink(final LineLink lineLink, boolean isInferred)
    {
        if(mLineData == null)
            mLineData = new LineData(LineLinker.hasPolyATail(this));

        if(isInferred)
            mLineData.InferredLinkedLineBreakends = lineLink;
        else
            mLineData.LinkedLineBreakends = lineLink;
    }

    public String toString()
    {
        return String.format("%s:%s pos(%s:%d:%d)", VcfId, type(), Chromosome, Position, Orient.asByte());
    }

    public String coordStr()
    {
        return String.format("%s:%d:%d", Context.getContig(), Position, Orient.asByte());
    }
}

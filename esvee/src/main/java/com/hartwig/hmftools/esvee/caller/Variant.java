package com.hartwig.hmftools.esvee.caller;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.sv.SvUtils.hasShortIndelLength;
import static com.hartwig.hmftools.common.sv.SvUtils.isIndel;
import static com.hartwig.hmftools.common.sv.SvVcfTags.ASM_INFO;
import static com.hartwig.hmftools.common.sv.SvVcfTags.AVG_FRAG_LENGTH;
import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.esvee.common.SvConstants.ASSEMBLY_INFO_DELIM;
import static com.hartwig.hmftools.esvee.common.SvConstants.JUNCTION_COORD_DELIM;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.esvee.assembly.types.Junction;
import com.hartwig.hmftools.esvee.caller.annotation.RepeatMaskAnnotation;
import com.hartwig.hmftools.esvee.common.FilterType;

import htsjdk.variant.variantcontext.VariantContext;

public class Variant
{
    // VCF data
    private final StructuralVariantType mType;

    private final Breakend[] mBreakends;

    // repeatedly used values for filtering are cached
    private final double mQual;
    private final String mInsertSequence;
    private boolean mIsShortLocalIndel;
    private int mPonCount;
    private boolean mIsHotspot;
    private boolean mGermline;
    private RepeatMaskAnnotation mRmAnnotation;
    private final Set<FilterType> mFilters;

    private List<Junction> mOriginalJunctions;

    public Variant(final StructuralVariant sv, final GenotypeIds genotypeIds)
    {
        if(isIndel(sv.type()))
        {
            mIsShortLocalIndel = hasShortIndelLength(abs(sv.end().position() - sv.start().position()));
        }
        else
        {
            mIsShortLocalIndel = false;
        }

        mType = sv.type();

        Breakend breakendStart = Breakend.from(
                this, true, sv.start(), sv.startContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);

        Breakend breakendEnd = sv.end() != null ?
                Breakend.from(this, false, sv.end(), sv.endContext(), genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal) : null;

        mBreakends = new Breakend[] { breakendStart, breakendEnd };

        mQual = sv.qualityScore();

        mInsertSequence = sv.insertSequence();

        mFilters = Sets.newHashSet();

        // keep any non-pass filters from assembly
        addExistingFilters(sv.startContext());
        addExistingFilters(sv.endContext());

        mPonCount = 0;
        mIsHotspot = false;
        mGermline = false;
        mRmAnnotation = null;
        mOriginalJunctions = null;
    }

    public String chromosomeStart() { return mBreakends[SE_START].Chromosome; }
    public String chromosomeEnd() { return !isSgl() ? mBreakends[SE_END].Chromosome : ""; }

    public int posStart() { return mBreakends[SE_START].Position; }
    public int posEnd() { return !isSgl() ? mBreakends[SE_END].Position : -1; }

    public Orientation orientStart() { return mBreakends[SE_START].Orient; }
    public Orientation orientEnd() { return !isSgl() ? mBreakends[SE_END].Orient : null; }

    public StructuralVariantType type() { return mType; }
    public boolean isSgl() { return mType == SGL; }

    public Breakend[] breakends() { return mBreakends; }
    public Breakend breakendStart() { return mBreakends[SE_START]; }
    public Breakend breakendEnd() { return mBreakends[SE_END]; }

    public VariantContext contextStart() { return mBreakends[SE_START].Context; }
    public VariantContext contextEnd() { return !isSgl() ? mBreakends[SE_END].Context : null; }

    public double qual() { return mQual; }
    public String insertSequence() { return mInsertSequence; }

    public int ponCount() { return mPonCount; }
    public void setPonCount(int count) { mPonCount = count; }

    public RepeatMaskAnnotation getRmAnnotation() { return mRmAnnotation; }
    public void setRepeatMaskAnnotation(final RepeatMaskAnnotation annotation) { mRmAnnotation = annotation; }

    public boolean isShortLocal() { return mIsShortLocalIndel; }

    public static boolean hasLength(final StructuralVariantType type)
    {
        return type == INS || type == INV || type == DEL || type == DUP;
    }

    public int svLength() { return hasLength(mType) ? abs(posEnd() - posStart()) : 0; }

    public int adjustedLength()
    {
        if(!hasLength(mType))
            return 0;

        if(mType == INS)
            return mInsertSequence.length();

        int positionLength = posEnd() - posStart();

        if(mType == DUP)
            ++positionLength;

        return positionLength + mInsertSequence.length();
    }

    public int averageFragmentLength() { return contextStart().getAttributeAsInt(AVG_FRAG_LENGTH, 0); }

    public int splitFragmentCount() { return contextStart().getAttributeAsInt(TOTAL_FRAGS, 0); }

    public void markHotspot() { mIsHotspot = true; }
    public boolean isHotspot() { return mIsHotspot; }

    public void markGermline() { mGermline = true; }
    public boolean isGermline() { return mGermline; }

    public boolean isLineSite()
    {
        for(Breakend breakend : mBreakends)
        {
            if(breakend != null && (breakend.isLine() || breakend.lineSiteBreakend() != null))
                return true;
        }

        return false;
    }

    public boolean inChainedAssembly()
    {
        for(Breakend breakend : mBreakends)
        {
            if(breakend != null && breakend.inChainedAssembly())
                return true;
        }

        return false;
    }

    public int maxUniqueFragmentPositions()
    {
        int maxUps = 0;
        for(Breakend breakend : mBreakends)
        {
            if(breakend != null)
                maxUps = max(maxUps, breakend.uniqueFragmentPositions());
        }

        return maxUps;
    }

    public void addFilter(final FilterType filter) { mFilters.add(filter); }
    public Set<FilterType> filters() { return mFilters; }
    public boolean isPass() { return mFilters.isEmpty(); }
    public boolean isFiltered() { return !mFilters.isEmpty(); }

    private void addExistingFilters(final VariantContext variantContext)
    {
        if(variantContext == null)
            return;

        // keep any non-pass filters from assembly
        for(String filterStr : variantContext.getFilters())
        {
            if(filterStr.equals(PASS_FILTER))
                continue;

            FilterType filterType = Arrays.stream(FilterType.values()).filter(x -> x.vcfTag().equals(x)).findFirst().orElse(null);

            if(filterType != null)
                mFilters.add(filterType);
        }
    }

    public List<Junction> originalJunctions()
    {
        if(mOriginalJunctions != null)
            return mOriginalJunctions;

        String assemblyInfoStr = mBreakends[0].Context.getAttributeAsString(ASM_INFO, "");

        if(assemblyInfoStr.isEmpty())
        {
            mOriginalJunctions = Collections.emptyList();
        }
        else
        {
            // chr10:87940542:1_chr10:87952584:-1
            String[] assemblyStrs = assemblyInfoStr.split(ASSEMBLY_INFO_DELIM, -1);
            mOriginalJunctions = Lists.newArrayListWithCapacity(assemblyStrs.length);

            for(String assemblyStr : assemblyStrs)
            {
                String[] items = assemblyStr.split(JUNCTION_COORD_DELIM);

                if(items.length < 3)
                    continue;

                String chromosome = items[0];
                int position = Integer.parseInt(items[1]);
                Orientation orientation = Orientation.fromByteStr(items[2]);

                mOriginalJunctions.add(new Junction(chromosome, position, orientation));
            }
        }

        return mOriginalJunctions;
    }

    public String toString()
    {
        if(!isSgl())
        {
            return String.format("%s pos(%s:%d - %s:%d)",
                    mType.toString(), chromosomeStart(), posStart(), chromosomeEnd(), posEnd());
        }
        else
        {
            return String.format("%s pos(%s:%d)",
                    mType.toString(), chromosomeStart(), posStart());
        }
    }
}

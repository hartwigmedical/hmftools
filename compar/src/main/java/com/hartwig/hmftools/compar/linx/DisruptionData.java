package com.hartwig.hmftools.compar.linx;

import static java.lang.String.format;

import static com.hartwig.hmftools.compar.common.Category.DISRUPTION;
import static com.hartwig.hmftools.compar.common.Category.values;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.MismatchType.VALUE;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.errorprone.annotations.Var;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;

public class DisruptionData implements ComparableItem
{
    public final String GeneName;

    public final List<StructuralVariantData> Variants;
    public final List<LinxBreakend> Breakends;

    protected static final String FLD_BREAKEND_INFO = "BreakendInfo";
    protected static final String FLD_REGION_TYPE = "RegionType";
    protected static final String FLD_CODING_CONTEXT = "CodingContext";
    protected static final String FLD_GENE_ORIENT = "GeneOrientation";
    protected static final String FLD_NEXT_SPLICE = "NextSpliceExonRank";

    public DisruptionData(
            final String geneName, final List<StructuralVariantData> variants, final List<LinxBreakend> breakends)
    {
        GeneName = geneName;
        Variants = variants;
        Breakends = breakends;
    }

    @Override
    public Category category() { return DISRUPTION; }

    @Override
    public String key()
    {
        return String.format("%s variants(%d) breakends(%d)",
                GeneName, Variants.size(), Breakends.size());
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();

        /*
        values.add(String.format("%s", Breakend.reportedDisruption()));
        values.add(String.format("%s", Breakend.regionType()));
        values.add(String.format("%s", Breakend.codingType()));
        values.add(String.format("%s", Breakend.geneOrientation()));
        values.add(String.format("%d", Breakend.nextSpliceExonRank()));
        */

        return values;
    }

    @Override
    public boolean reportable() { return Breakends.stream().anyMatch(x -> x.reportedDisruption()); }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final DisruptionData otherDisruptionData = (DisruptionData)other;

        return GeneName.equals(otherDisruptionData.GeneName);

        /*
        if(otherSv.SvData.type() != SvData.type())
            return false;

        if(!otherSv.SvData.startChromosome().equals(mComparisonPositionStart.Chromosome)
        || !otherSv.SvData.endChromosome().equals(mComparisonPositionEnd.Chromosome))
            return false;

        if(otherSv.SvData.startPosition() != mComparisonPositionStart.Position
        || otherSv.SvData.endPosition() != mComparisonPositionEnd.Position)
            return false;

        if(otherSv.SvData.startOrientation() != SvData.startOrientation() || otherSv.SvData.endOrientation() != SvData.endOrientation())
            return false;

        if(!otherSv.Breakend.gene().equals(Breakend.gene()))
            return false;

        if((otherSv.mCheckTranscript || mCheckTranscript) && !otherSv.Breakend.transcriptId().equals(Breakend.transcriptId()))
            return false;

        return true;
        */
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds)
    {
        final DisruptionData otherDisruptionData = (DisruptionData)other;

        final List<String> diffs = Lists.newArrayList();

        // compare each breakend and record differences
        List<LinxBreakend> breakends = Lists.newArrayList(Breakends);
        List<LinxBreakend> otherBreakends = Lists.newArrayList(otherDisruptionData.Breakends);

        int index = 0;
        while(index < breakends.size())
        {
            LinxBreakend breakend = breakends.get(index);
            StructuralVariantData variant = findVariant(breakend);

            LinxBreakend otherBreakend = findMatchingBreakend(variant, breakend);

            if(otherBreakend != null)
            {
                checkDiff(diffs, FLD_REGION_TYPE, breakend.regionType().toString(), otherBreakend.regionType().toString());
                checkDiff(diffs, FLD_CODING_CONTEXT, breakend.codingType().toString(), otherBreakend.codingType().toString());
                checkDiff(diffs, FLD_REPORTED, breakend.reportedDisruption(), otherBreakend.reportedDisruption());
                checkDiff(diffs, FLD_NEXT_SPLICE, breakend.nextSpliceExonRank(), otherBreakend.nextSpliceExonRank());

                breakends.remove(index);
                otherBreakends.remove(otherBreakend);
            }
            else
            {
                // record an unmatched breakend or SV
                diffs.add(format("unmatched SV(%s %s:%d-%s-%d)",
                        variant.type(), variant.startChromosome(), variant.endChromosome(), variant.endPosition()));

                ++index;
            }
        }

        return !diffs.isEmpty() ? new Mismatch(this, other, VALUE, diffs) : null;
    }

    public StructuralVariantData findVariant(final LinxBreakend breakend)
    {
        return Variants.stream().filter(x -> x.id() == breakend.svId()).findFirst().orElse(null);
    }

    public LinxBreakend findMatchingBreakend(final StructuralVariantData otherVariant, final LinxBreakend otherBreakend)
    {
        for(StructuralVariantData variant : Variants)
        {
            if(otherVariant.type() != variant.type())
                continue;

            if(!otherVariant.startChromosome().equals(variant.startChromosome()) || !otherVariant.endChromosome().equals(variant.endChromosome()))
                continue;

            // could match within homology
            if(otherVariant.startPosition() != otherVariant.startPosition() || otherVariant.endPosition() != variant.endPosition())
                continue;

            if(otherVariant.startOrientation() != variant.startOrientation() || otherVariant.endOrientation() != variant.endOrientation())
                continue;

            // now search breakends
            List<LinxBreakend> breakends = Breakends.stream().filter(x -> x.svId() == variant.id()).collect(Collectors.toList());

            for(LinxBreakend breakend : breakends)
            {
                if(breakend.transcriptId().equals(otherBreakend.transcriptId()))
                    return breakend;
            }
        }

        return null;
    }
}

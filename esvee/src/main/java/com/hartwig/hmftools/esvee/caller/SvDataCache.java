package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.esvee.caller.FilterConstants.PANEL_INCLUSION_BUFFER;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.GenotypeIds;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;

public class SvDataCache
{
    private final CallerConfig mConfig;
    private final TargetRegions mTargetRegions;
    private final StructuralVariantFactory mSvFactory;

    private final List<Variant> mSvData;
    private final Map<String,List<Breakend>> mChromosomeBreakends;

    private int mHardFilteredCount;

    public SvDataCache(final CallerConfig config, final TargetRegions targetRegions)
    {
        mConfig = config;
        mSvData = Lists.newArrayList();
        mChromosomeBreakends = Maps.newHashMap();

        mTargetRegions = targetRegions;

        mSvFactory = new StructuralVariantFactory(new CompoundFilter(false));
        mHardFilteredCount = 0;
    }

    public void setGenotypeOrdinals(final GenotypeIds genotypeIds)
    {
        mSvFactory.setGenotypeOrdinals(genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);
    }

    public List<Variant> getSvList() { return mSvData; }
    public Map<String,List<Breakend>> getBreakendMap() { return mChromosomeBreakends; }

    public int sglCount() { return (int)mSvData.stream().filter(x -> x.isSgl()).count(); }
    public int svCount() { return (int)mSvData.stream().filter(x -> !x.isSgl()).count(); }
    public int hardFilteredCount() { return mHardFilteredCount; }
    public int incompleteSVs() { return mSvFactory.unmatched().size(); }

    public void processVariant(final VariantContext variant, final GenotypeIds genotypeIds)
    {
        if(!HumanChromosome.contains(variant.getContig()))
            return;

        boolean isSgl = StructuralVariantFactory.isSingleBreakend(variant);

        if(isSgl)
        {
            if(mTargetRegions.hasTargetRegions() && !mTargetRegions.inTargetRegions(variant.getContig(), variant.getStart()))
            {
                ++mHardFilteredCount;
                return;
            }

            StructuralVariant sv = mSvFactory.createSingleBreakend(variant);
            addSvData(new Variant(sv, genotypeIds));
            return;
        }

        String mateId = StructuralVariantFactory.mateId(variant);

        if(mateId == null)
            return;

        int currentSvCount = mSvFactory.results().size();
        mSvFactory.addVariantContext(variant);

        // check if both breakends have now been encountered
        if(currentSvCount == mSvFactory.results().size())
            return;

        StructuralVariant sv = popLastSv(); // get and clear from storage

        if(sv == null)
            return;

        // one of the breakends at least must be within a targeted region
        if(mTargetRegions.hasTargetRegions() && !svWithinTargetedRegion(sv))
            return;

        addSvData(new Variant(sv, genotypeIds));
    }

    private boolean svWithinTargetedRegion(final StructuralVariant sv)
    {
        return mTargetRegions.inTargetRegions(sv.chromosome(true), sv.position(true), PANEL_INCLUSION_BUFFER)
            || mTargetRegions.inTargetRegions(sv.chromosome(false), sv.position(false), PANEL_INCLUSION_BUFFER);
    }

    private StructuralVariant popLastSv()
    {
        if(mSvFactory.results().isEmpty())
            return null;

        StructuralVariant sv = mSvFactory.results().get(0);
        mSvFactory.results().remove(0);

        return sv;
    }

    private void addSvData(final Variant var)
    {
        // optionally filter out by config
        if(mConfig.excludeVariant(var))
            return;

        mSvData.add(var);
    }

    public void buildBreakendMap()
    {
        buildBreakendMap(mSvData, mChromosomeBreakends);
    }

    public static void buildBreakendMap(final List<Variant> variants, final Map<String,List<Breakend>> chrBreakendMap)
    {
        for(Variant var : variants)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                Breakend breakend = var.breakends()[se];

                if(breakend == null)
                    continue;

                List<Breakend> breakends = chrBreakendMap.get(breakend.Chromosome);

                if(breakends == null)
                {
                    breakends = Lists.newArrayList();
                    chrBreakendMap.put(breakend.Chromosome, breakends);
                }

                breakends.add(breakend);
            }
        }

        for(List<Breakend> breakends : chrBreakendMap.values())
        {
            Collections.sort(breakends, new BreakendPositionComparator());
        }
    }

    public static class BreakendPositionComparator implements Comparator<Breakend>
    {
        public int compare(final Breakend first, final Breakend second)
        {
            if(first.Position == second.Position)
            {
                if(first.Orient == second.Orient)
                    return 0;
                else
                    return first.Orient.isForward() ? -1 : 1;
            }
            else
            {
                return first.Position < second.Position ? -1 : 1;
            }
        }
    }

    public void clear()
    {
        mChromosomeBreakends.clear();
        mSvData.clear();
    }
}

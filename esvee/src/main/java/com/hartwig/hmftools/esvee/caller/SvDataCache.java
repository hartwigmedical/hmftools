package com.hartwig.hmftools.esvee.caller;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

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

        final StructuralVariant sv = popLastSv(); // get and clear from storage

        if(sv == null)
            return;

        addSvData(new Variant(sv, genotypeIds));
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
        for(Variant var : mSvData)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                Breakend breakend = var.breakends()[se];

                if(breakend == null)
                    continue;

                List<Breakend> breakends = mChromosomeBreakends.get(breakend.Chromosome);

                if(breakends == null)
                {
                    breakends = Lists.newArrayList();
                    mChromosomeBreakends.put(breakend.Chromosome, breakends);
                }

                int index = 0;
                while(index < breakends.size())
                {
                    if(breakend.Position < breakends.get(index).Position)
                        break;

                    ++index;
                }

                breakends.add(index, breakend);
            }
        }

        for(List<Breakend> breakends : mChromosomeBreakends.values())
        {
            for(int index = 0; index < breakends.size(); ++index)
            {
                breakends.get(index).setChrLocationIndex(index);
            }
        }
    }

    public List<Breakend> selectOthersNearby(final Breakend breakend, int additionalDistance, int maxSeekDistance)
    {
        List<Breakend> breakends = mChromosomeBreakends.get(breakend.Chromosome);

        List<Breakend> closeBreakends = Lists.newArrayList();

        if(breakends == null)
            return closeBreakends;

        int minStart = breakend.minPosition() - additionalDistance;
        int maxStart = breakend.maxPosition() + additionalDistance;

        // search down
        for(int index = breakend.chrLocationIndex() - 1; index >= 0; --index)
        {
            Breakend nextBreakend = breakends.get(index);

            if(nextBreakend.sv() == breakend.sv())
                continue;

            if(nextBreakend.maxPosition() < breakend.minPosition() - maxSeekDistance)
                break;

            if(nextBreakend.minPosition() <= maxStart && nextBreakend.maxPosition() >= minStart)
                closeBreakends.add(0, nextBreakend);
        }

        for(int index = breakend.chrLocationIndex() + 1; index < breakends.size(); ++index)
        {
            Breakend nextBreakend = breakends.get(index);

            if(nextBreakend.sv() == breakend.sv())
                continue;

            if(nextBreakend.minPosition() > breakend.maxPosition() + maxSeekDistance)
                break;

            if(nextBreakend.minPosition() <= maxStart && nextBreakend.maxPosition() >= minStart)
                closeBreakends.add(nextBreakend);
        }

        return closeBreakends;
    }

    public void clear()
    {
        mChromosomeBreakends.clear();
        mSvData.clear();
    }
}

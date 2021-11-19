package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.SvData;

import org.apache.commons.compress.utils.Lists;

public class SvDataCache
{
    private final List<SvData> mSvData;
    private final Map<String,List<Breakend>> mChromosomeBreakends;

    public SvDataCache()
    {
        mSvData = Lists.newArrayList();
        mChromosomeBreakends = Maps.newHashMap();
    }

    public List<SvData> getSvList() { return mSvData; }

    public void addSvData(final SvData sv) { mSvData.add(sv); }

    public void buildBreakendMap()
    {
        for(SvData sv : mSvData)
        {
            // keep SGL breakends?
            for(int se = SE_START; se <= SE_END; ++se)
            {
                Breakend breakend = sv.breakends()[se];

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
}

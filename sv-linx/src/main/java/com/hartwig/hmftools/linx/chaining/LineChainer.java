package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.NO_LINE_ELEMENT;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.types.LinkType.TEMPLATED_INSERTION;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.compress.utils.Lists;

public class LineChainer
{
    private int mClusterId;
    private final List<SvVarData> mSvList;
    private final List<SvLinkedPair> mAssembledLinks;

    private final List<SvChain> mChains;
    private final List<SvBreakend> mSourceBreakends;
    private final Set<String> mSourceChromosomes;

    public LineChainer()
    {
        mClusterId = -1;
        mSvList = Lists.newArrayList();
        mAssembledLinks = Lists.newArrayList();
        mChains = Lists.newArrayList();
        mSourceBreakends = Lists.newArrayList();
        mSourceChromosomes = Sets.newHashSet();
    }

    public final List<SvChain> getChains() { return mChains; }

    public void clear()
    {
        mClusterId = -1;
        mSvList.clear();
        mChains.clear();
        mAssembledLinks.clear();
        mSourceBreakends.clear();
        mSourceChromosomes.clear();
    }

    public void initialise(SvCluster cluster)
    {
        clear();
        mClusterId = cluster.id();
        mSvList.addAll(cluster.getSVs());
        mAssembledLinks.addAll(cluster.getAssemblyLinkedPairs());
    }

    public void formChains()
    {
        // in many cases the source elements will have already formed into assembled links to the insertion site
        // and these should be removed from consideration and stored in independent chains
        for(SvVarData var : mSvList)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                final SvBreakend breakend = var.getBreakend(se);

                if(breakend == null || var.getLineElement(isStart(se)).equals(NO_LINE_ELEMENT))
                    continue;

                // skip breakends already assembled into links
                if(mAssembledLinks.stream().anyMatch(x -> x.hasBreakend(breakend)))
                    continue;

                mSourceBreakends.add(breakend);
                mSourceChromosomes.add(breakend.chromosome());
            }
        }

        // form chains from existing assembled links
        for(SvLinkedPair pair : mAssembledLinks)
        {
            if(tryAddLink(pair))
                continue;

            addNewChain(pair);
        }

        if(!mSourceBreakends.isEmpty())
        {
            addSourceBreakends();
        }
    }

    private void addNewChain(final SvLinkedPair pair)
    {
        SvChain newChain = new SvChain(mChains.size());
        newChain.addLink(pair, true);
        mChains.add(newChain);

        LNX_LOGGER.debug("new chain({}) with pair({})", newChain.id(), pair);
    }

    private boolean tryAddLink(final SvLinkedPair pair)
    {
        for(SvChain chain : mChains)
        {
            // only add links on source chromosomes, not at insertion sites
            for(int se = SE_START; se <= SE_END; ++se)
            {
                SvBreakend chainBreakend = chain.getOpenBreakend(se);

                if(chainBreakend == null || !mSourceChromosomes.contains(chainBreakend.chromosome()))
                    continue;

                if(chain.canAddLinkedPair(pair, isStart(se), false))
                {
                    LNX_LOGGER.debug("adding pair({}) to existing chain({}))", pair, chain.id());
                    chain.addLink(pair, isStart(se));
                    return true;
                }
            }
        }

        return false;
    }

    private void addSourceBreakends()
    {
        while(true)
        {
            final List<SvLinkedPair> possiblePairs = Lists.newArrayList();
            SvLinkedPair shortestPair = null;

            for(int i = 0; i < mSourceBreakends.size(); ++i)
            {
                SvBreakend breakend = mSourceBreakends.get(i);

                // form links amongst the source breakends
                for(int j = i + 1; j < mSourceBreakends.size(); ++j)
                {
                    SvBreakend breakend2 = mSourceBreakends.get(j);

                    SvLinkedPair newPair = tryFormLinkedPair(breakend, breakend2);

                    if(newPair == null)
                        continue;

                    possiblePairs.add(newPair);

                    if(shortestPair == null || newPair.length() < shortestPair.length())
                        shortestPair = newPair;
                }

                for(SvChain chain : mChains)
                {
                    // only add links on source chromosomes, not at insertion sites
                    for(int se = SE_START; se <= SE_END; ++se)
                    {
                        SvBreakend chainBreakend = chain.getOpenBreakend(se);

                        if(chainBreakend == null || !mSourceChromosomes.contains(chainBreakend.chromosome()))
                            continue;

                        SvLinkedPair newPair = tryFormLinkedPair(breakend, chainBreakend);

                        if(newPair == null)
                            continue;

                        possiblePairs.add(newPair);

                        if(shortestPair == null || newPair.length() < shortestPair.length())
                            shortestPair = newPair;
                    }
                }
            }

            if(shortestPair == null)
                return;

            if(!tryAddLink(shortestPair))
            {
                addNewChain(shortestPair);
            }

            mSourceBreakends.remove(shortestPair.firstBreakend());
            mSourceBreakends.remove(shortestPair.secondBreakend());

            if(mSourceBreakends.isEmpty())
                return;
        }
    }

    private static SvLinkedPair tryFormLinkedPair(final SvBreakend breakend1, final SvBreakend breakend2)
    {
        if(breakend1 == null || breakend2 == null)
            return null;

        if(!LinkFinder.areLinkedSection(breakend1.getSV(), breakend2.getSV(), breakend1.usesStart(), breakend2.usesStart()))
            return null;

        int minTILength = getMinTemplatedInsertionLength(breakend2, breakend1);
        if(abs(breakend2.position() - breakend1.position()) < minTILength)
            return null;

        return new SvLinkedPair(
                breakend1.getSV(), breakend2.getSV(), TEMPLATED_INSERTION, breakend1.usesStart(), breakend2.usesStart());
    }

}

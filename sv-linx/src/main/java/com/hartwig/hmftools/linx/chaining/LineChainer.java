package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.isPossibleLink;
import static com.hartwig.hmftools.linx.types.LinkType.TEMPLATED_INSERTION;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_TEMPLATED_INSERTION_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.compress.utils.Lists;

public class LineChainer
{
    private int mClusterId;
    private String mSampleId;
    private final List<SvVarData> mSvList;
    private final List<SvLinkedPair> mAssembledLinks;

    private final List<SvChain> mChains;
    private final List<SvBreakend> mSourceBreakends;
    private final Set<String> mSourceChromosomes;

    private BufferedWriter mFileWriter;

    public LineChainer()
    {
        mClusterId = -1;
        mSampleId = "";
        mSvList = Lists.newArrayList();
        mAssembledLinks = Lists.newArrayList();
        mChains = Lists.newArrayList();
        mSourceBreakends = Lists.newArrayList();
        mSourceChromosomes = Sets.newHashSet();
        mFileWriter = null;
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

    public void initialise(final String sampleId, final SvCluster cluster)
    {
        clear();
        mClusterId = cluster.id();
        mSampleId = sampleId;
        mSvList.addAll(cluster.getSVs());
        mAssembledLinks.addAll(cluster.getAssemblyLinkedPairs());
    }

    public void formChains()
    {
        // in many cases the source elements will have already formed into assembled links to the insertion site
        // and these should be removed from consideration and stored in independent chains

        final List<SvVarData> sglMapped = Lists.newArrayList();

        for(SvVarData var : mSvList)
        {
            for(int se = SE_START; se <= SE_END; ++se)
            {
                final SvBreakend breakend = var.getBreakend(se);

                if(breakend == null)
                {
                    //  create local breakends if the alignment mappings suggest a link could be made
                    if(!var.getSglMappings().isEmpty())
                        sglMapped.add(var);

                    continue;
                }

                if(!var.isLineElement(isStart(se)))
                    continue;

                mSourceChromosomes.add(breakend.chromosome());

                // skip breakends already assembled into links
                if(mAssembledLinks.stream().anyMatch(x -> x.hasBreakend(breakend)))
                    continue;

                mSourceBreakends.add(breakend);
            }
        }

        addSglMappings(sglMapped);

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

        writeChainData();
        writeUnchainedSVs();
    }

    private void addNewChain(final SvLinkedPair pair)
    {
        SvChain newChain = new SvChain(mChains.size());

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final SvBreakend breakend = pair.getBreakend(se);
            if(breakend.getSV().isSglBreakend() && !breakend.usesStart())
            {
                breakend.getSV().setSglEndBreakend(breakend);
            }
        }

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
        // put any remaining breakends into chains by making inferred TIs, but only if either the new linked pair has the other
        // breakends going to insertion sites, or the new linked pair joins 2 chains with their other ends going similarly
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

    private SvLinkedPair tryFormLinkedPair(final SvBreakend breakend1, final SvBreakend breakend2)
    {
        if(breakend1 == null || breakend2 == null)
            return null;

        int minTILength = getMinTemplatedInsertionLength(breakend2, breakend1);

        if(!isPossibleLink(
                breakend1.chromosome(), breakend1.position(), breakend1.orientation(),
                breakend2.chromosome(), breakend2.position(), breakend2.orientation(), minTILength))
        {
            return null;
        }

        // the other ends of these breakends must either go to the same insertion location or be part of chains which do the same
        final SvBreakend otherBreakend1 = getOtherBreakend(breakend1);
        final SvBreakend otherBreakend2 = getOtherBreakend(breakend2);

        if(otherBreakend1 == null || otherBreakend2 == null)
            return null;

        if(otherBreakend1.getDBLink() == null || otherBreakend1.getDBLink() != otherBreakend2.getDBLink())
            return null;

        return new SvLinkedPair(breakend1, breakend2, TEMPLATED_INSERTION);
    }

    private SvBreakend getOtherBreakend(final SvBreakend breakend)
    {
        for(SvChain chain : mChains)
        {
            // only add links on source chromosomes, not at insertion sites
            for(int se = SE_START; se <= SE_END; ++se)
            {
                SvBreakend chainBreakend = chain.getOpenBreakend(se);

                if(chainBreakend == breakend)
                    return chain.getOpenBreakend(switchIndex(se));
            }
        }

        return breakend.getOtherBreakend();
    }

    private void addSglMappings(final List<SvVarData> sglMapped)
    {
        final Set<SvVarData> sglCandidates = Sets.newHashSet();

        if(mSourceBreakends.isEmpty())
        {
            for(int i = 0; i < sglMapped.size(); ++i)
            {
                final SvVarData sgl1 = sglMapped.get(i);

                if(sglCandidates.contains(sgl1))
                    continue;

                for(int j = i+1; j < sglMapped.size(); ++j)
                {
                    final SvVarData sgl2 = sglMapped.get(j);

                    if(sglCandidates.contains(sgl2))
                        continue;

                    for(final SglMapping mapping1 : sgl1.getSglMappings())
                    {
                        final SglMapping mapping2 = sgl2.getSglMappings().stream()
                                .filter(x -> x.possibleLink(mapping1)).findFirst().orElse(null);

                        if(mapping2 != null)
                        {
                            mSourceBreakends.add(new SvBreakend(sgl1, mapping1));
                            mSourceBreakends.add(new SvBreakend(sgl2, mapping2));
                            sglCandidates.add(sgl1);
                            sglCandidates.add(sgl2);
                            break;
                        }
                    }
                }
            }
        }
        else
        {
            for(SvVarData sgl : sglMapped)
            {
                boolean linked = false;

                for(SvBreakend breakend : mSourceBreakends)
                {
                    int minDistance = max(breakend.anchorDistance() - breakend.homology().length(), MIN_TEMPLATED_INSERTION_LENGTH);

                    for(final SglMapping mapping : sgl.getSglMappings())
                    {
                        if(isPossibleLink(mapping.Chromosome, mapping.Position, mapping.Orientation,
                                breakend.chromosome(), breakend.position(), breakend.orientation(), minDistance))
                        {
                            mSourceBreakends.add(new SvBreakend(sgl, mapping));
                            linked = true;
                            break;
                        }
                    }

                    if(linked)
                        break;
                }
            }
        }
    }

    public void initialiseOutput(final String outputDir)
    {
        if(outputDir == null)
            return;

        try
        {
            String outputFileName = outputDir + "LNX_LINE_CHAINS.csv";
            mFileWriter = createBufferedWriter(outputFileName, false);

            mFileWriter.write("SampleId,ClusterId,ChainId,ChainSvCount,AsmbLinks,ChainDesc,SourceChr,SourcePosStart,SourcePosEnd");
            mFileWriter.write(",InsertChr,InsertPosStart,InsertPosEnd,SourceInvPosStart,SourceInvPosEnd");
            mFileWriter.newLine();
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error initialising line chain data: {}", e.toString());
        }
    }

    private void writeChainData()
    {
        if(mFileWriter == null || mChains.isEmpty())
            return;

        try
        {
            for(final SvChain chain : mChains)
            {
                final SvBreakend chainStart = chain.getOpenBreakend(true);
                final SvBreakend chainEnd = chain.getOpenBreakend(false);

                final int[] typeCounts = new int[StructuralVariantType.values().length];

                for(final SvVarData var : chain.getSvList())
                {
                    ++typeCounts[typeAsInt(var.type())];
                }

                mFileWriter.write(String.format("%s,%d,%d,%d,%d,%s",
                        mSampleId, mClusterId, chain.id(), chain.getSvCount(), chain.getAssemblyLinkCount(), getSvTypesStr(typeCounts)));

                final SvLinkedPair firstLink = chain.getLinkedPairs().get(0);
                final SvLinkedPair lastLink = chain.getLinkedPairs().get(chain.getLinkedPairs().size() - 1);

                mFileWriter.write(String.format(",%s,%d,%d",
                        firstLink.chromosome(), firstLink.firstBreakend().position(), lastLink.secondBreakend().position()));

                mFileWriter.write(String.format(",%s,%d,%d",
                        chainStart != null ? chainStart.chromosome() : "-1",
                        chainStart != null ? chainStart.position() : -1,
                        chainEnd != null ? chainEnd.position() : -1));

                final SvVarData inv = chain.getSvList().stream().filter(x -> x.type() == INV).findFirst().orElse(null);

                if(inv != null)
                {
                    mFileWriter.write(String.format(",%d,%d", inv.position(true), inv.position(false)));
                }
                else
                {
                    mFileWriter.write(",-1,-1");
                }

                mFileWriter.newLine();

            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing line chain data: {}", e.toString());
        }
    }

    private void writeUnchainedSVs()
    {
        if(mFileWriter == null || mSourceBreakends.isEmpty())
            return;

        try
        {
            int nonChainId = mChains.size();

            for(final SvBreakend breakend : mSourceBreakends)
            {
                final SvBreakend otherBreakend = breakend.getOtherBreakend();
                boolean isLineInv = breakend.getSV().type() == INV && breakend.getSV().inLineElement();

                if(mSourceBreakends.contains(otherBreakend))
                {
                    if(breakend.position() > otherBreakend.position())
                        continue;
                }

                mFileWriter.write(String.format("%s,%d,%d,%d,%d,%s,%s",
                        mSampleId, mClusterId, nonChainId++, 1, 0, breakend.getSV().type(), breakend.chromosome()));

                if(breakend.inLineElement() && !isLineInv)
                {
                    mFileWriter.write(String.format(",%d,%d",
                            breakend.position(), otherBreakend != null && otherBreakend.inLineElement() ? otherBreakend.position() : -1));
                }
                else
                {
                    mFileWriter.write(",-1,-1");
                }

                if(otherBreakend != null && !otherBreakend.inLineElement() && !isLineInv)
                {
                    mFileWriter.write(String.format(",%s,%d,-1",
                            otherBreakend.chromosome(), otherBreakend.position()));
                }
                else
                {
                    mFileWriter.write(",-1,-1,-1");
                }

                if(isLineInv)
                {
                    mFileWriter.write(String.format(",%d,%d", breakend.position(), otherBreakend.position()));
                }
                else
                {
                    mFileWriter.write(",-1,-1");
                }

                mFileWriter.newLine();

            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing line chain data: {}", e.toString());
        }
    }

    public void close() { closeBufferedWriter(mFileWriter); }

}

package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.switchIndex;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.typeAsInt;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getSvTypesStr;
import static com.hartwig.hmftools.linx.annotators.LineElementAnnotator.LINE_ELEMENT_PROXIMITY_DISTANCE;
import static com.hartwig.hmftools.linx.annotators.LineElementType.SUSPECT;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.isPossibleLink;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_DEL_LENGTH;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_TEMPLATED_INSERTION_LENGTH;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.CohortFileInterface;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;

import org.apache.commons.compress.utils.Lists;

public class LineChainer implements CohortFileInterface
{
    private int mClusterId;
    private String mSampleId;
    private final CohortDataWriter mCohortDataWriter;
    private boolean mLoggingEnabled;
    private final List<SvVarData> mSvList;
    private final List<LinkedPair> mAssembledLinks;

    private final List<SvChain> mChains;
    private final List<SvBreakend> mSourceBreakends;
    private final Set<String> mSourceChromosomes;

    private final Set<SvVarData> mLoggedSVs;

    public LineChainer(final CohortDataWriter cohortDataWriter)
    {
        mLoggingEnabled = false;
        mCohortDataWriter = cohortDataWriter;
        mClusterId = -1;
        mSampleId = "";
        mSvList = Lists.newArrayList();
        mAssembledLinks = Lists.newArrayList();
        mChains = Lists.newArrayList();
        mSourceBreakends = Lists.newArrayList();
        mLoggedSVs = Sets.newHashSet();
        mSourceChromosomes = Sets.newHashSet();
    }

    public final List<SvChain> getChains() { return mChains; }

    public void enableLogging() { mLoggingEnabled = true; }

    public void clear()
    {
        mClusterId = -1;
        mSvList.clear();
        mChains.clear();
        mAssembledLinks.clear();
        mSourceBreakends.clear();
        mSourceChromosomes.clear();
        mLoggedSVs.clear();
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
        for(LinkedPair pair : mAssembledLinks)
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

    private void addNewChain(final LinkedPair pair)
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

    private boolean tryAddLink(final LinkedPair pair)
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
            final List<LinkedPair> possiblePairs = Lists.newArrayList();
            LinkedPair shortestPair = null;

            for(int i = 0; i < mSourceBreakends.size(); ++i)
            {
                SvBreakend breakend = mSourceBreakends.get(i);

                // form links amongst the source breakends
                for(int j = i + 1; j < mSourceBreakends.size(); ++j)
                {
                    SvBreakend breakend2 = mSourceBreakends.get(j);

                    LinkedPair newPair = tryFormLinkedPair(breakend, breakend2);

                    if(newPair == null)
                        continue;

                    possiblePairs.add(newPair);

                    if(shortestPair == null || newPair.positionDistance() < shortestPair.positionDistance())
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

                        LinkedPair newPair = tryFormLinkedPair(breakend, chainBreakend);

                        if(newPair == null)
                            continue;

                        possiblePairs.add(newPair);

                        if(shortestPair == null || newPair.positionDistance() < shortestPair.positionDistance())
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

    private LinkedPair tryFormLinkedPair(final SvBreakend breakend1, final SvBreakend breakend2)
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

        // if too far away are likely in different line elements
        if(abs(breakend1.position() - breakend2.position()) > LINE_ELEMENT_PROXIMITY_DISTANCE * 2)
            return null;

        // the other ends of these breakends must either go to the same insertion location or be part of chains which do the same
        final SvBreakend otherBreakend1 = getOtherBreakend(breakend1);
        final SvBreakend otherBreakend2 = getOtherBreakend(breakend2);

        if(otherBreakend1 == null || otherBreakend2 == null)
            return null;

        if(otherBreakend1.getDBLink() == null || otherBreakend1.getDBLink() != otherBreakend2.getDBLink())
            return null;

        return new LinkedPair(breakend1, breakend2);
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

                            // mark as line elements so they'll feature in the chain data logged
                            sgl1.addLineElement(SUSPECT, false);
                            sgl2.addLineElement(SUSPECT, false);
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
                            sgl.addLineElement(SUSPECT, false);
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

    private static final String COHORT_WRITER_LINE_CHAINER = "LineChainer";

    @Override
    public String fileType() { return COHORT_WRITER_LINE_CHAINER; }

    @Override
    public BufferedWriter createWriter(final String outputDir)
    {
        try
        {
            String outputFileName = outputDir + "LNX_LINE_CHAINS.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            // definitional fields
            writer.write("SampleId,ClusterId,ChainId,ChainSvCount,AsmbLinks,ChainDesc,ChainComplete");
            writer.write(",SourceChr,SourcePosStart,SourcePosEnd,SourceOrientStart,SourceOrientEnd");
            writer.write(",InsertChr,InsertPosStart,InsertPosEnd,SourceInvPosStart,SourceInvPosEnd,SourceInvOrient");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to initialise line chainer: {}", e.toString());
            return null;
        }
    }

    private void writeChainData()
    {
        if(mChains.isEmpty() || !mLoggingEnabled || mCohortDataWriter == null)
            return;

        List<String> outputLines = Lists.newArrayList();

        for(final SvChain chain : mChains)
        {
            final int[] typeCounts = new int[StructuralVariantType.values().length];

            for(final SvVarData var : chain.getSvList())
            {
                ++typeCounts[typeAsInt(var.type())];
            }

            final SvBreakend chainStart = chain.getOpenBreakend(true);
            final SvBreakend chainEnd = chain.getOpenBreakend(false);
            final LinkedPair firstLink = chain.getLinkedPairs().get(0);
            final LinkedPair lastLink = chain.getLinkedPairs().get(chain.getLinkedPairs().size() - 1);

            int chainSvCount = chain.getSvCount();
            final int[] sourcePositions = { -1, -1 };
            final byte[] sourceOrients = { 0, 0 };
            int lowerSourceIndex = 0;
            int lowerInsertIndex = 0;
            int[] insertPositions = { -1, -1 };
            String insertChr = "0";
            final int[] invPositions = { -1, -1 };
            byte invOrient = 0;

            chain.getSvList().stream().forEach(x -> mLoggedSVs.add(x));

            final SvVarData inv = chain.getSvList().stream()
                    .filter(x -> x.type() == INV)
                    .filter(x -> x.inLineElement())
                    .findFirst().orElse(null);

            if(firstLink.firstBreakend().inLineElement() && firstLink.firstBreakend().type() != INV)
            {
                sourcePositions[0] = firstLink.firstBreakend().position();
                sourceOrients[0] = firstLink.firstBreakend().orientation();
            }

            if(lastLink.secondBreakend().inLineElement() && lastLink.secondBreakend().type() != INV)
            {
                if(sourcePositions[0] == -1)
                {
                    sourcePositions[0] = lastLink.secondBreakend().position();
                    sourceOrients[0] = lastLink.secondBreakend().orientation();
                }
                else
                {
                    sourcePositions[1] = lastLink.secondBreakend().position();
                    sourceOrients[1] = lastLink.secondBreakend().orientation();
                }
            }

            if(sourcePositions[0] > 0 && sourcePositions[1] > 0)
                lowerSourceIndex = sourcePositions[0] <= sourcePositions[1] ? 0 : 1;

            if(chainStart != null && !chainStart.inLineElement())
            {
                insertChr = chainStart.chromosome();
                insertPositions[0] = chainStart.position();
            }

            if(chainEnd != null && !chainEnd.inLineElement())
            {
                insertChr = chainEnd.chromosome();
                insertPositions[1] = chainEnd.position();
            }

            // check for an unmapped SGL or INF in a DB at the insertion site
            if(insertPositions[1] < 0 && insertPositions[0] > 0 && chainStart.getDBLink() != null)
            {
                final DbPair dbPair = chainStart.getDBLink();
                final SvBreakend otherBreakend = dbPair.getOtherBreakend(chainStart);
                final SvVarData otherSV = otherBreakend.getSV();
                if(dbPair.length() <= MIN_DEL_LENGTH &&
                        mSvList.contains(otherSV) && !chain.getSvList().contains(otherSV))
                {
                    insertPositions[1] = otherBreakend.position();
                    mLoggedSVs.add(otherSV);
                    ++typeCounts[typeAsInt(otherSV.type())];
                    ++chainSvCount;
                }
            }
            else if(insertPositions[0] < 0 && insertPositions[1] > 0 && chainEnd.getDBLink() != null)
            {
                final DbPair dbPair = chainEnd.getDBLink();
                final SvBreakend otherBreakend = dbPair.getOtherBreakend(chainEnd);
                final SvVarData otherSV = otherBreakend.getSV();
                if(dbPair.length() <= MIN_DEL_LENGTH && mSvList.contains(otherSV) && !chain.getSvList().contains(otherSV))
                {
                    insertPositions[0] = otherBreakend.position();
                    mLoggedSVs.add(otherSV);
                    ++typeCounts[typeAsInt(otherSV.type())];
                    ++chainSvCount;
                }
            }

            if(insertPositions[0] > 0 && insertPositions[1] > 0)
                lowerInsertIndex = insertPositions[0] <= insertPositions[1] ? 0 : 1;

            if(inv != null)
            {
                invPositions[SE_START] = inv.position(true);
                invPositions[SE_END] = inv.position(false);
                invOrient = inv.orientation(true);
            }

            StringBuilder sb = new StringBuilder();

            sb.append(String.format("%s,%d,%d,%d,%d,%s,true",
                    mSampleId, mClusterId, chain.id(), chainSvCount, chain.getAssemblyLinkCount(), getSvTypesStr(typeCounts)));

            sb.append(String.format(",%s,%d,%d,%d,%d",
                    firstLink.chromosome(), sourcePositions[lowerSourceIndex], sourcePositions[switchIndex(lowerSourceIndex)],
                    sourceOrients[lowerSourceIndex], sourceOrients[switchIndex(lowerSourceIndex)]));

            sb.append(String.format(",%s,%d,%d",
                    insertChr, insertPositions[lowerInsertIndex], insertPositions[switchIndex(lowerInsertIndex)]));

            sb.append(String.format(",%d,%d,%d",
                    invPositions[SE_START], invPositions[SE_END], invOrient));

            outputLines.add(sb.toString());
        }

        mCohortDataWriter.write(this, outputLines);
    }

    private void writeUnchainedSVs()
    {
        if(!mLoggingEnabled || mCohortDataWriter == null)
            return;

        final List<SvVarData> unchainedSVs = mSvList.stream().filter(x -> !mLoggedSVs.contains(x)).collect(Collectors.toList());

        int nonChainId = mChains.size();

        List<String> outputLines = Lists.newArrayList();

        for(final SvVarData var : unchainedSVs)
        {
            if(mLoggedSVs.contains(var))
                continue;

            final int[] typeCounts = new int[StructuralVariantType.values().length];
            ++typeCounts[typeAsInt(var.type())];

            String sourceChr = "";
            final int[] sourcePositions = { -1, -1 };
            int lowerSourceIndex = 0;
            int lowerInsertIndex = 0;
            final int[] insertPositions = { -1, -1 };
            final int[] invPositions = { -1, -1 };
            byte invOrient = 0;
            String insertChr = "0";
            byte[] sourceOrient = {0, 0};
            int chainSvCount = 1;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                final SvBreakend breakend = var.getBreakend(se);
                if(breakend == null)
                    continue;

                if(breakend.inLineElement())
                {
                    if(breakend.type() == INV && breakend.getOtherBreakend().inLineElement())
                    {
                        invOrient = breakend.orientation();
                        invPositions[se] = breakend.position();
                    }
                    else
                    {
                        sourceChr = breakend.chromosome();
                        sourcePositions[se] = breakend.position();
                        sourceOrient[se] = breakend.orientation();
                    }
                }
                else
                {
                    insertChr = breakend.chromosome();
                    insertPositions[se] = breakend.position();

                    // factor in 2 breakends in a remote DB but not chained at their source site
                    // typically SGL=2, BND=1_SGL=1, BND=2 (diff orients) or BND=1_INV=1_SGL=1
                    if(breakend.getDBLink() != null)
                    {
                        final DbPair dbPair = breakend.getDBLink();
                        final SvBreakend otherBreakend = dbPair.getOtherBreakend(breakend);
                        final SvVarData otherSV = otherBreakend.getSV();
                        if(dbPair.length() <= MIN_DEL_LENGTH && mSvList.contains(otherSV) && !mLoggedSVs.contains(otherSV))
                        {
                            // check for a matching source line element
                            final SvBreakend firstOtherBreakend = breakend.getOtherBreakend();
                            final SvBreakend secondOtherBreakend = otherBreakend.getOtherBreakend();

                            boolean isValidPairing = (firstOtherBreakend == null || secondOtherBreakend == null)
                                    || (firstOtherBreakend.inLineElement() && secondOtherBreakend.inLineElement()
                                    && firstOtherBreakend.chromosome().equals(secondOtherBreakend.chromosome())
                                    && abs(firstOtherBreakend.position() - secondOtherBreakend.position()) < LINE_ELEMENT_PROXIMITY_DISTANCE);

                            if(isValidPairing)
                            {
                                insertPositions[switchIndex(se)] = otherBreakend.position();
                                mLoggedSVs.add(otherSV);
                                ++typeCounts[typeAsInt(otherSV.type())];
                                ++chainSvCount;

                                if(firstOtherBreakend != null)
                                {
                                    sourceChr = firstOtherBreakend.chromosome();
                                    sourcePositions[0] = firstOtherBreakend.position();
                                    sourceOrient[0] = firstOtherBreakend.orientation();
                                }

                                if(secondOtherBreakend != null)
                                {
                                    sourceChr = secondOtherBreakend.chromosome();
                                    sourcePositions[1] = secondOtherBreakend.position();
                                    sourceOrient[1] = secondOtherBreakend.orientation();
                                }

                                break;
                            }
                        }
                    }
                }
            }

            if(sourcePositions[0] > 0 && sourcePositions[1] > 0)
                lowerSourceIndex = sourcePositions[0] <= sourcePositions[1] ? 0 : 1;

            if(insertPositions[0] > 0 && insertPositions[1] > 0)
                lowerInsertIndex = insertPositions[0] <= insertPositions[1] ? 0 : 1;

            // SampleId,ClusterId,ChainId,ChainSvCount,AsmbLinks,ChainDesc,SourceChr,SourcePosStart,SourcePosEnd");
            // InsertChr,InsertPosStart,InsertPosEnd,SourceInvPosStart,SourceInvPosEnd,SourceInvOrient,LoneBeOrient");
            StringBuilder sb = new StringBuilder();

            sb.append(String.format("%s,%d,%d,%d,%d,%s,false",
                    mSampleId, mClusterId, nonChainId++, chainSvCount, 0, getSvTypesStr(typeCounts)));

            sb.append(String.format(",%s,%d,%d,%d,%d,%s,%d,%d",
                    sourceChr, sourcePositions[lowerSourceIndex], sourcePositions[switchIndex(lowerSourceIndex)],
                    sourceOrient[lowerSourceIndex], sourceOrient[switchIndex(lowerSourceIndex)],
                    insertChr, insertPositions[lowerInsertIndex], insertPositions[switchIndex(lowerInsertIndex)]));

            sb.append(String.format(",%d,%d,%d",
                    invPositions[SE_START], invPositions[SE_END], invOrient));

            outputLines.add(sb.toString());
        }

        mCohortDataWriter.write(this, outputLines);
    }

}

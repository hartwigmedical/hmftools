package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INS;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.types.SvConstants.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_DB;
import static com.hartwig.hmftools.linx.types.SvLinkedPair.LINK_TYPE_TI;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.types.SvLinkedPair;

public class LinkFinder
{
    private boolean mLogVerbose;

    public LinkFinder()
    {
        mLogVerbose = false;
    }

    public void setLogVerbose(boolean toggle)
    {
        mLogVerbose = toggle;
    }

    public void findAssembledLinks(SvCluster cluster)
    {
        cluster.setAssemblyLinkedPairs(createAssemblyLinkedPairs(cluster));
    }

    public static List<SvLinkedPair> createAssemblyLinkedPairs(SvCluster cluster)
    {
        List<SvLinkedPair> linkedPairs = Lists.newArrayList();

        // find 2 breakends with matching assembly info and form them into a linked pair
        if(cluster.getSvCount() < 2)
            return linkedPairs;

        boolean hasMultiBreakendLinks = false;
        final Map<String,List<SvBreakend>> chrBreakendMap = cluster.getChrBreakendMap();

        List<SvBreakend> multiConnectionBreakends = Lists.newArrayList();

        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            if(breakendList.size() == 1)
                continue;

            for (int i = 0; i < breakendList.size() -1; ++i)
            {
                final SvBreakend lowerBreakend = breakendList.get(i);

                if(lowerBreakend.orientation() != -1)
                    continue;

                final SvVarData lowerSV = lowerBreakend.getSV();

                if(lowerSV.type() == INS || lowerSV.isSglBreakend())
                    continue;

                for (int j = i+1; j < breakendList.size(); ++j)
                {
                    final SvBreakend upperBreakend = breakendList.get(j);

                    if(upperBreakend.orientation() != 1)
                        continue;

                    final SvVarData upperSV = upperBreakend.getSV();

                    if(upperBreakend.getSV() == lowerBreakend.getSV())
                        continue;

                    if(upperSV.type() == INS || upperSV.isSglBreakend())
                        continue;

                    boolean v1Start = lowerBreakend.usesStart();
                    boolean v2Start = upperBreakend.usesStart();

                    if (!haveLinkedAssemblies(lowerSV, upperSV, v1Start, v2Start))
                        continue;

                    // it's possible for a breakend to already be in an assembled link
                    // eg the replicated breakend scenario in COLO829's BND from chr 3-6
                    // but allow this if its ploidy supports this

                    if(lowerBreakend.isAssembledLink() && !multiConnectionBreakends.contains(lowerBreakend))
                    {
                        multiConnectionBreakends.add(lowerBreakend);
                        hasMultiBreakendLinks = true;
                    }

                    if(upperBreakend.isAssembledLink() && !multiConnectionBreakends.contains(upperBreakend))
                    {
                        multiConnectionBreakends.add(upperBreakend);
                        hasMultiBreakendLinks = true;
                    }

                    // form a new TI from these 2 BEs
                    SvLinkedPair newPair = new SvLinkedPair(lowerSV, upperSV, LINK_TYPE_TI, v1Start, v2Start);
                    newPair.setIsAssembled();
                    lowerSV.addLinkedPair(newPair, v1Start);
                    upperSV.addLinkedPair(newPair, v2Start);

                    linkedPairs.add(newPair);

                    // to avoid logging unlikely long TIs
                    LNX_LOGGER.debug("cluster({}) adding assembly linked {} pair({}) length({})",
                            cluster.id(), newPair.linkType(), newPair.toString(), newPair.length());
                }
            }
        }

        if(hasMultiBreakendLinks)
        {
            hasMultiBreakendLinks = analyseMultiConnectionBreakends(cluster, multiConnectionBreakends, linkedPairs);

            if(hasMultiBreakendLinks)
                cluster.setRequiresReplication();
        }

        return linkedPairs;
    }

    private static boolean analyseMultiConnectionBreakends(
            final SvCluster cluster, final List<SvBreakend> multiConnectionBreakends, List<SvLinkedPair> linkedPairs)
    {
        // check for the scenario where one link A-B is made but also A-C and B-C, where B is a short simple SV
        boolean reassessMultiConnection = false;

        for(final SvBreakend breakend : multiConnectionBreakends)
        {
            final List<SvLinkedPair> links = breakend.getSV().getAssembledLinkedPairs(breakend.usesStart());

            if(links.size() < 2)
                continue;

            List<SvLinkedPair> spanningLinks = Lists.newArrayList();

            for(final SvLinkedPair pair : links)
            {
                if(spanningLinks.contains(pair))
                    continue;

                final SvBreakend otherBreakend = pair.getOtherBreakend(breakend);
                final SvVarData otherSv = otherBreakend.getSV();

                if(otherSv.type() != DEL && otherSv.type() != DUP)
                    continue;

                // now look for a common breakend partner
                final List<SvLinkedPair> otherLinks = otherSv.getAssembledLinkedPairs(!otherBreakend.usesStart());

                if(otherLinks.isEmpty())
                    continue;

                for(final SvLinkedPair otherPair : links)
                {
                    if(otherPair == pair)
                        continue;

                    final SvBreakend otherBreakend1 = otherPair.getOtherBreakend(breakend);

                    if(otherLinks.stream().anyMatch(x -> x.hasBreakend(otherBreakend1)))
                    {
                        LNX_LOGGER.debug("cluster({}) breakends({} & {} type={}) both assemble to breakend({}), removing outer link({})",
                                cluster.id(), breakend, otherBreakend, otherSv.type(), otherBreakend1, otherPair);

                        spanningLinks.add(otherPair);
                        reassessMultiConnection = true;
                    }
                }
            }

            for(SvLinkedPair pair : spanningLinks)
            {
                breakend.getSV().getLinkedPairs(breakend.usesStart()).remove(pair);

                final SvBreakend otherBreakend = pair.getOtherBreakend(breakend);
                otherBreakend.getSV().getLinkedPairs(otherBreakend.usesStart()).remove(pair);

                linkedPairs.remove(pair);
            }
        }

        boolean hasMultiConnection = true;

        if(reassessMultiConnection)
        {
            hasMultiConnection = false;

            for(final SvLinkedPair pair : linkedPairs)
            {
                for(int se = SE_START; se <= SE_END; ++se)
                {
                    final SvBreakend breakend = pair.getBreakend(se);

                    if(linkedPairs.stream().filter(x -> x != pair).anyMatch(x -> x.hasBreakend(breakend)))
                    {
                        hasMultiConnection = true;
                        break;
                    }
                }
            }
        }

        return hasMultiConnection;
    }

    public static boolean haveLinkedAssemblies(final SvVarData var1, final SvVarData var2, boolean v1Start, boolean v2Start)
    {
        return var1.getTIAssemblies(v1Start).stream().anyMatch(x -> var2.getTIAssemblies(v2Start).contains(x));
    }

    public static boolean areLinkedSection(final SvVarData v1, final SvVarData v2, boolean v1Start, boolean v2Start)
    {
        // templated insertions are allowed to traverse the centromere
        if(v1.position(v1Start) < 0 || v2.position(v2Start) < 0)
            return false;

        if(!v1.chromosome(v1Start).equals(v2.chromosome(v2Start)))
            return false;

        // start apart and heading towards each other
        long pos1 = v1.position(v1Start);
        boolean headsLeft1 = (v1.orientation(v1Start) == 1);
        long pos2 = v2.position(v2Start);
        boolean headsLeft2 = (v2.orientation(v2Start) == 1);

        boolean breakendsFace = false;
        if(pos1 < pos2 && !headsLeft1 && headsLeft2)
            breakendsFace = true;
        else if(pos2 < pos1 && headsLeft1 && !headsLeft2)
            breakendsFace = true;

        if(!breakendsFace)
            return false;

        return true;
    }

    public static int getMinTemplatedInsertionLength(SvBreakend breakend1, SvBreakend breakend2)
    {
        // return the minimum length from each breakend's anchor distances if present
        int homologyLengths = breakend1.homology().length() + breakend2.homology().length();
        int minAnchorDistance = min(breakend1.anchorDistance(), breakend2.anchorDistance());
        return max(minAnchorDistance - homologyLengths, MIN_TEMPLATED_INSERTION_LENGTH);
    }

    public static void findDeletionBridges(final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size() - 1; ++i)
            {
                SvBreakend breakend = breakendList.get(i);
                SvBreakend nextBreakend = breakendList.get(i+1);
                SvVarData var1 = breakend.getSV();
                SvVarData var2 = nextBreakend.getSV();

                if(var1 == var2)
                    continue;

                if(breakend.orientation() == nextBreakend.orientation())
                    continue;

                if(breakend.arm() != nextBreakend.arm())
                    continue;

                long distance = nextBreakend.position() - breakend.position();
                int minTiLength = getMinTemplatedInsertionLength(breakend, nextBreakend);

                if(breakend.orientation() == 1 && nextBreakend.orientation() == -1)
                {
                    // breakends face away as per a normal DB
                    SvLinkedPair dbPair = new SvLinkedPair(var1, var2, LINK_TYPE_DB, breakend.usesStart(), nextBreakend.usesStart());
                    markDeletionBridge(dbPair);
                }
                else if(distance < minTiLength)
                {
                    // facing breakends with a distance less than the anchor distances are in fact a DB with overlap
                    if(haveLinkedAssemblies(var1, var2, breakend.usesStart(), nextBreakend.usesStart()))
                    {
                        // however check that they don't have assemblies on these breakends with each other
                        continue;
                    }

                    SvLinkedPair dbPair = new SvLinkedPair(var1, var2, LINK_TYPE_TI, breakend.usesStart(), nextBreakend.usesStart());
                    markDeletionBridge(dbPair);
                }
            }
        }
    }

    private static void markDeletionBridge(SvLinkedPair dbPair)
    {
        // if either SV has a shorter DB already in existence, then keep it
        final SvLinkedPair existingDBLink1 = dbPair.first().getDBLink(dbPair.firstLinkOnStart());
        final SvLinkedPair existingDBLink2 = dbPair.second().getDBLink(dbPair.secondLinkOnStart());

        if((existingDBLink1 != null && existingDBLink1.length() < dbPair.length())
        || (existingDBLink2 != null && existingDBLink2.length() < dbPair.length()))
        {
            return;
        }

        // remove any conflicting DB link info
        if(existingDBLink1 != null)
        {
            if(existingDBLink1.first() == dbPair.first())
                existingDBLink1.second().setDBLink(null, existingDBLink1.secondLinkOnStart());
            else
                existingDBLink1.first().setDBLink(null, existingDBLink1.firstLinkOnStart());
        }

        if(existingDBLink2 != null)
        {
            if(existingDBLink2.first() == dbPair.first())
                existingDBLink2.second().setDBLink(null, existingDBLink2.secondLinkOnStart());
            else
                existingDBLink2.first().setDBLink(null, existingDBLink2.firstLinkOnStart());
        }

        dbPair.first().setDBLink(dbPair, dbPair.firstLinkOnStart());
        dbPair.second().setDBLink(dbPair, dbPair.secondLinkOnStart());
    }

    public List<SvLinkedPair> createInferredLinkedPairs(SvCluster cluster, List<SvVarData> svList, boolean allowSingleBEs)
    {
        List<SvLinkedPair> linkedPairs = Lists.newArrayList();

        for (int i = 0; i < svList.size(); ++i)
        {
            SvVarData var1 = svList.get(i);

            if(var1.type() == INS || (var1.isSglBreakend() && !allowSingleBEs))
                continue;

            for(int be1 = SE_START; be1 <= SE_END; ++be1)
            {
                boolean v1Start = isStart(be1);

                if(var1.isSglBreakend() && !v1Start)
                    continue;

                // if an assembly linked pair has already been created for this breakend, look no further
                if(var1.hasAssemblyLink(v1Start))
                    continue;

                for (int j = i+1; j < svList.size(); ++j)
                {
                    SvVarData var2 = svList.get(j);

                    if(var1 == var2)
                        continue;

                    if(var2.type() == INS || (var2.isSglBreakend() && !allowSingleBEs))
                        continue;

                    for(int be2 = SE_START; be2 <= SE_END; ++be2)
                    {
                        boolean v2Start = isStart(be2);

                        if(var1.isSglBreakend() && !v1Start)
                            continue;

                        if(var2.hasAssemblyLink(v2Start))
                            continue;

                        long distance = abs(var1.position(v1Start) - var2.position(v2Start));
                        int minTiLength = getMinTemplatedInsertionLength(var1.getBreakend(v1Start), var2.getBreakend(v2Start));

                        if (!areLinkedSection(var1, var2, v1Start, v2Start) || distance < minTiLength)
                            continue;

                        // form a new TI from these 2 BEs
                        SvLinkedPair newPair = new SvLinkedPair(var1, var2, LINK_TYPE_TI, v1Start, v2Start);

                        if(newPair.linkType() == LINK_TYPE_DB)
                        {
                            // was considered too short to be a TI and so converted to a DB
                            continue;
                        }

                        // insert in order of increasing length
                        int index = 0;
                        boolean skipNewPair = false;

                        for (; index < linkedPairs.size(); ++index)
                        {
                            final SvLinkedPair pair = linkedPairs.get(index);

                            // check for matching BEs on a pair that is much shorter, and if so skip creating this new linked pair
                            if(newPair.length() >= pair.length() && newPair.hasLinkClash(pair))
                            {
                                // allow a TI if only a DB has been found
                                skipNewPair = true;
                                break;
                            }

                            if (pair.length() > newPair.length())
                                break;
                        }

                        if(skipNewPair)
                            continue;

                        if (index >= linkedPairs.size())
                            linkedPairs.add(newPair);
                        else
                            linkedPairs.add(index, newPair);

                        // if(newPair.length() < mUtils.getBaseDistance())
                        if(mLogVerbose) //  && !var1.isReplicatedSv() && !var2.isReplicatedSv()
                        {
                            // to avoid logging unlikely long TIs
                            LNX_LOGGER.debug("cluster({}) adding inferred linked {} pair({}) length({}) at index({})",
                                    cluster.id(), newPair.linkType(), newPair.toString(), newPair.length(), index);
                        }
                    }
                }
            }
        }

        // prior to consolidating linked pairs, check for duplicate BE in the spanning SVs
        // matchDuplicateBEToLinkedPairs(linkedPairs, spanningSVs);

        if(linkedPairs.isEmpty())
            return linkedPairs;

        LNX_LOGGER.debug("cluster({}) has {} inferred linked pairs", cluster.id(), linkedPairs.size());

        return linkedPairs;
    }
}

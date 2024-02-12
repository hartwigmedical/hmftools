package com.hartwig.hmftools.linx.chaining;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INF;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.types.LinxConstants.MIN_TEMPLATED_INSERTION_LENGTH;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.types.LinkedPair;

public class LinkFinder
{
    public static List<LinkedPair> createAssemblyLinkedPairs(SvCluster cluster)
    {
        List<LinkedPair> linkedPairs = Lists.newArrayList();

        // find 2 breakends with matching assembly info and form them into a linked pair
        if(cluster.getSvCount() < 2)
            return linkedPairs;

        boolean hasMultiBreakendLinks = false;
        final Map<String,List<SvBreakend>> chrBreakendMap = cluster.getChrBreakendMap();

        List<SvBreakend> multiConnectionBreakends = Lists.newArrayList();

        for(final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            if(breakendList.size() == 1)
                continue;

            for(int i = 0; i < breakendList.size() -1; ++i)
            {
                final SvBreakend lowerBreakend = breakendList.get(i);

                if(lowerBreakend.orientation() != -1)
                    continue;

                final SvVarData lowerSV = lowerBreakend.getSV();

                if(lowerSV.type() == INS || lowerSV.isSglBreakend())
                    continue;

                for(int j = i+1; j < breakendList.size(); ++j)
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

                    if(!haveLinkedAssemblies(lowerSV, upperSV, v1Start, v2Start))
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
                    LinkedPair newPair = new LinkedPair(lowerBreakend, upperBreakend);
                    newPair.setIsAssembled();
                    lowerSV.addLinkedPair(newPair, v1Start);
                    upperSV.addLinkedPair(newPair, v2Start);

                    linkedPairs.add(newPair);

                    LNX_LOGGER.debug("cluster({}) adding assembly linked pair({}) length({})",
                            cluster.id(), newPair.toString(), newPair.baseLength());
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
            final SvCluster cluster, final List<SvBreakend> multiConnectionBreakends, List<LinkedPair> linkedPairs)
    {
        // check for the scenario where one link A-B is made but also A-C and B-C, where B is a short simple SV
        boolean reassessMultiConnection = false;

        for(final SvBreakend breakend : multiConnectionBreakends)
        {
            final List<LinkedPair> links = breakend.getSV().getAssembledLinkedPairs(breakend.usesStart());

            if(links.size() < 2)
                continue;

            List<LinkedPair> spanningLinks = Lists.newArrayList();

            for(final LinkedPair pair : links)
            {
                if(spanningLinks.contains(pair))
                    continue;

                final SvBreakend otherBreakend = pair.getOtherBreakend(breakend);
                final SvVarData otherSv = otherBreakend.getSV();

                if(otherSv.type() != DEL && otherSv.type() != DUP)
                    continue;

                // now look for a common breakend partner
                final List<LinkedPair> otherLinks = otherSv.getAssembledLinkedPairs(!otherBreakend.usesStart());

                if(otherLinks.isEmpty())
                    continue;

                for(final LinkedPair otherPair : links)
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

            for(LinkedPair pair : spanningLinks)
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

            for(final LinkedPair pair : linkedPairs)
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

    public static boolean isPossibleLink(final String chr1, int pos1, byte orient1, final String chr2, int pos2, byte orient2, int minDistance)
    {
        if(!chr1.equals(chr2))
            return false;

        return isPossibleLink(pos1, orient1, pos2, orient2, minDistance);
    }

    public static boolean isPossibleLink(int pos1, byte orient1, int pos2, byte orient2, int minDistance)
    {
        if(orient1 == orient2)
            return false;

        return (pos1 <= pos2 - minDistance && orient1 == NEG_ORIENT) || (pos1 >= pos2 + minDistance && orient1 == POS_ORIENT);
    }

    public static int getMinTemplatedInsertionLength(SvBreakend breakend1, SvBreakend breakend2)
    {
        if(breakend1.type() == INF || breakend2.type() == INF)
        {
            if(abs(breakend1.position() - breakend2.position()) <= 1)
                return 0;
        }

        // return the minimum length from each breakend's anchor distances if present
        int homologyLengths = breakend1.homology().length() + breakend2.homology().length();
        int minAnchorDistance = min(breakend1.anchorDistance(), breakend2.anchorDistance());
        return max(minAnchorDistance - homologyLengths, MIN_TEMPLATED_INSERTION_LENGTH);
    }

    public static void findDeletionBridges(final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        for(final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            for(int i = 0; i < breakendList.size() - 1; ++i)
            {
                SvBreakend breakend = breakendList.get(i);
                SvBreakend nextBreakend = breakendList.get(i + 1);
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
                    // skip if the next breakend is a short overlapping DB
                    if(i < breakendList.size() - 2)
                    {
                        SvBreakend nextNextBreakend = breakendList.get(i + 2);
                        int nextDistance = nextNextBreakend.position() - nextBreakend.position();
                        int nextMinTiLength = getMinTemplatedInsertionLength(nextBreakend, nextNextBreakend);

                        if(nextDistance < distance && (nextDistance < nextMinTiLength || nextMinTiLength == 0))
                            continue;
                    }

                    // breakends face away as per a normal DB
                    DbPair dbPair = new DbPair(breakend, nextBreakend);
                    markDeletionBridge(dbPair);
                }
                else if(distance < minTiLength || minTiLength == 0) // factoring in INFs at the same base
                {
                    // facing breakends with a distance less than the anchor distances are in fact a DB with overlap
                    if(haveLinkedAssemblies(var1, var2, breakend.usesStart(), nextBreakend.usesStart()))
                    {
                        // however check that they don't have assemblies on these breakends with each other
                        continue;
                    }

                    DbPair dbPair = new DbPair(breakend, nextBreakend);
                    markDeletionBridge(dbPair);
                }
            }
        }
    }

    private static void markDeletionBridge(DbPair dbPair)
    {
        // if either SV has a shorter DB already in existence, then keep it
        final DbPair existingDBLink1 = dbPair.lowerSV().getDBLink(dbPair.lowerLinkOnStart());
        final DbPair existingDBLink2 = dbPair.upperSV().getDBLink(dbPair.upperLinkOnStart());

        if((existingDBLink1 != null && existingDBLink1.length() < dbPair.length())
        || (existingDBLink2 != null && existingDBLink2.length() < dbPair.length()))
        {
            return;
        }

        // remove any conflicting DB link info
        if(existingDBLink1 != null)
        {
            if(existingDBLink1.lowerSV() == dbPair.lowerSV())
                existingDBLink1.upperSV().setDBLink(null, existingDBLink1.upperLinkOnStart());
            else
                existingDBLink1.lowerSV().setDBLink(null, existingDBLink1.lowerLinkOnStart());
        }

        if(existingDBLink2 != null)
        {
            if(existingDBLink2.lowerSV() == dbPair.lowerSV())
                existingDBLink2.upperSV().setDBLink(null, existingDBLink2.upperLinkOnStart());
            else
                existingDBLink2.lowerSV().setDBLink(null, existingDBLink2.lowerLinkOnStart());
        }

        dbPair.lowerSV().setDBLink(dbPair, dbPair.lowerLinkOnStart());
        dbPair.upperSV().setDBLink(dbPair, dbPair.upperLinkOnStart());
    }

}

package com.hartwig.hmftools.linx.analysis;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.INV;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.CHAIN_ASSEMBLY_LINK_COUNT;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.CHAIN_LENGTH;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.CHAIN_LINK_COUNT;
import static com.hartwig.hmftools.linx.chaining.ChainUtils.breakendsAreChained;
import static com.hartwig.hmftools.linx.chaining.LinkFinder.getMinTemplatedInsertionLength;
import static com.hartwig.hmftools.linx.types.LinxConstants.MAX_FOLDBACK_CHAIN_LENGTH;

import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class FoldbackFinder
{
    public static void markFoldbacks(final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        markFoldbacks(chrBreakendMap, false);
    }

    public static boolean markFoldbacks(final Map<String, List<SvBreakend>> chrBreakendMap, boolean recheckChainedSVs)
    {
        boolean foundAnyFoldback = false;

        // find all valid consecutive breakends formed either from a single SV or a chained set
        for(final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            List<SvBreakend> breakendList = entry.getValue();

            for(int i = 0; i < breakendList.size() - 1; ++i)
            {
                SvBreakend breakend = breakendList.get(i);

                if(breakend.isAssembledLink())
                    continue;

                if(recheckChainedSVs && (breakend.isFoldback() || breakend.getCluster().getSvCount() == 1))
                    continue;

                SvBreakend nextBreakend = null;

                SvBreakend beFront = null; // the lower position for orientation +1 and vice versa
                SvBreakend beBack = null;

                int j = i + 1;
                while(j < breakendList.size())
                {
                    nextBreakend = breakendList.get(j);

                    // first skip over any breakends in a DB with the initial breakend
                    if(j == i + 1 && breakend.orientation() == -1 && nextBreakend.orientation() == 1
                    && nextBreakend.position() - breakend.position() < getMinTemplatedInsertionLength(nextBreakend, breakend))
                    {
                        ++j;
                        continue;
                    }

                    // check for any assembled links in between the potential foldback breakends
                    if(j + 1 < breakendList.size() && nextBreakend.isAssembledLink()
                    && nextBreakend.getSV().getLinkedPair(nextBreakend.usesStart()) != null)
                    {
                        LinkedPair asmbLink = nextBreakend.getSV().getLinkedPair(nextBreakend.usesStart());
                        SvBreakend nextNextBreakend = breakendList.get(j + 1);
                        if(asmbLink.getOtherBreakend(nextBreakend) == nextNextBreakend)
                        {
                            // skip over both these assembled links
                            j += 2;
                            continue;
                        }
                    }

                    // check again for an overlapping DB at the outer (potential) foldback breakend
                    if(breakend.orientation() == 1 && nextBreakend.orientation() == -1 && j < breakendList.size() - 1)
                    {
                        SvBreakend nextNextBreakend = breakendList.get(j + 1);

                        if(nextNextBreakend.orientation() == breakend.orientation()
                        && nextNextBreakend.position() - nextBreakend.position() < getMinTemplatedInsertionLength(nextBreakend, nextNextBreakend))
                        {
                            nextBreakend = nextNextBreakend;
                        }
                    }

                    // now check for opposite orientation for the potential foldback
                    if(nextBreakend.orientation() == breakend.orientation() && !nextBreakend.isAssembledLink())
                    {
                        beFront = breakend.orientation() == 1 ? breakend : nextBreakend;
                        beBack = breakend.orientation() == 1 ? nextBreakend : breakend;
                    }

                    break;
                }

                boolean foldbackFound = false;

                if(beFront != null && beBack != null)
                {
                    // the foldback is invalid if it has a deletion bridge with overhang on the front-facing breakend
                    final DbPair dbLink = beFront.getSV().getDBLink(beFront.usesStart());

                    if(dbLink == null || dbLink.length() > 0)
                    {
                        if(!recheckChainedSVs || (beFront.getCluster() == beBack.getCluster() && !beFront.isFoldback() && !beBack.isFoldback()))
                        {
                            foldbackFound = checkFoldbackBreakends(beFront, beBack);
                        }
                    }
                }

                if(!foldbackFound && !recheckChainedSVs)
                {
                    checkReplicatedBreakendFoldback(breakend);
                }

                if(foldbackFound)
                    foundAnyFoldback = true;
            }
        }

        return foundAnyFoldback;
    }

    private static boolean checkFoldbackBreakends(SvBreakend beStart, SvBreakend beEnd)
    {
        // SVs are marked as being in a foldback if they are consecutive breakends,
        // have the same orientation, and are either an INV or part of a chain

        // beStart is the one with the lower position

        SvVarData varEnd = beEnd.getSV();
        SvVarData varStart = beStart.getSV();

        if(varEnd.type() == INS || varStart.type() == INS)
            return false;

        final SvCluster cluster = varEnd.getCluster();

        boolean singleSV = varEnd.equals(varStart);

        if(singleSV)
        {
            if(varStart.type() != INV)
                return false;
        }
        else
        {
            // must be same cluster
            if(varStart.getCluster() != cluster)
                return false;
        }

        // skip unclustered DELs & DUPs, reciprocal INV or reciprocal BNDs
        if(cluster.isResolved())
            return false;

        boolean beEndUsesStart = beEnd.usesStart();
        boolean beStartUsesStart = beStart.usesStart();

        String chainInfo = "0;0;0";

        if(singleSV)
        {
            // constraint is that the ends of this INV don't link to BND taking the path off this chromosome
            final SvChain chain = cluster.findChain(varEnd);

            if(chain != null)
            {
                int bndLinks = 0;
                for (final LinkedPair pair : chain.getLinkedPairs())
                {
                    if (pair.first() == varEnd && pair.second().type() == BND)
                        ++bndLinks;
                    else if (pair.second() == varEnd && pair.first().type() == BND)
                        ++bndLinks;

                    if (bndLinks == 2)
                        return false;
                }
            }
        }
        else
        {
            final SvChain chain1 = cluster.findChain(varEnd);
            final SvChain chain2 = cluster.findChain(varStart);

            if(chain1 == null || chain2 == null || chain1 != chain2)
                return false;

            // check if a path can be walked between these 2 breakends along the chain
            // without going back through this foldback point
            int[] chainData = breakendsAreChained(chain1, varEnd, !beEndUsesStart, varStart, !beStartUsesStart);

            if(chainData[CHAIN_LINK_COUNT] == 0 ) // || chainData[CHAIN_LINK_COUNT] != chainData[CHAIN_ASSEMBLY_LINK_COUNT]
                return false;

            int chainLength = chainData[CHAIN_LENGTH];

            if(chainLength > MAX_FOLDBACK_CHAIN_LENGTH)
            {
                /*
                LNX_LOGGER.info("sample({}) chained foldback breakends({} and {}) have long length({} links={} asmb={})",
                        mSampleId, beEnd.toString(), beStart.toString(),
                        chainLength, chainData[CHAIN_LINK_COUNT], chainData[CHAIN_ASSEMBLY_LINK_COUNT]);
                */
                return false;
            }

            chainInfo = String.format("%d;%d;%d",
                    chainData[CHAIN_LINK_COUNT], chainData[CHAIN_ASSEMBLY_LINK_COUNT], chainLength);
        }

        int length = abs(beEnd.position() - beStart.position());

        // if either variant already has foldback info set, favour
        // a) simple inversions then
        // b) shortest length

        boolean skipFoldback = false;
        if(varEnd.getFoldbackBreakend(beEndUsesStart) != null && !singleSV && varEnd.getFoldbackLength(beEndUsesStart) < length)
        {
            skipFoldback = true;
        }
        else if(varStart.getFoldbackBreakend(beStartUsesStart) != null && !singleSV && varStart.getFoldbackLength(beStartUsesStart) < length)
        {
            skipFoldback = true;
        }

        if(skipFoldback)
            return false;

        if(varEnd.getFoldbackBreakend(beEndUsesStart) != null)
        {
            final SvBreakend otherBreakend = varEnd.getFoldbackBreakend(beEndUsesStart);
            otherBreakend.getSV().setFoldbackLink(otherBreakend.usesStart(), null, -1, "");
        }

        if(varStart.getFoldbackBreakend(beStartUsesStart) != null)
        {
            final SvBreakend otherBreakend = varStart.getFoldbackBreakend(beStartUsesStart);
            otherBreakend.getSV().setFoldbackLink(otherBreakend.usesStart(), null, -1, "");
        }

        varEnd.setFoldbackLink(beEndUsesStart, beStart, length, chainInfo);
        varStart.setFoldbackLink(beStartUsesStart, beEnd, length, chainInfo);

        if(varEnd.equals(varStart))
        {
            LNX_LOGGER.debug("cluster({}) foldback inversion SV({}) length({})", cluster.id(), varEnd.posId(), length);
        }
        else
        {
            LNX_LOGGER.debug("cluster({}) foldback be1({}) be2({}) length({})", cluster.id(), beEnd.toString(), beStart.toString(), length);
        }

        return true;
    }

    private static void checkReplicatedBreakendFoldback(SvBreakend be)
    {
        // a special case where one ends of an SV connects to both ends of a single other variant
        // during a replication event and in doing so forms a foldback
        final SvVarData var = be.getSV();

        // if the other breakend is assembled to 2 other breakends and forms a chain with this breakend
        // open at both ends, then consider it a single-SV foldback
        if(var.getAssembledLinkedPairs(!be.usesStart()).size() < 2)
            return;

        final SvChain chain = var.getCluster().findChain(var);

        if(chain == null)
            return;

        if(chain.getOpenBreakend(true) == be && chain.getOpenBreakend(false) == be)
        {
            final String chainInfo = String.format("%d;%d;%d", 1, 1, chain.getLength(false));

            var.setFoldbackLink(be.usesStart(), be, 0, chainInfo);

            LNX_LOGGER.debug("cluster({}) foldback SV({} : {}) with own breakend({}) chain({})",
                    var.getCluster().id(), var.posId(), var.type(), be.toString(), chain.id());
        }
    }
}

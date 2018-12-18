package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.svanalysis.analysis.SvClusteringMethods.DEFAULT_PROXIMITY_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.breakendsMatch;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.getProximity;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.isSpecificCluster;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.haveLinkedAssemblies;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isSpecificSV;
import static com.hartwig.hmftools.svanalysis.types.SvVarData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_DIFF;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_DB;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class LinkFinder
{
    private boolean mLogVerbose;

    public static int MIN_TEMPLATED_INSERTION_LENGTH = 30;
    public static int NO_DB_MARKER = -(MIN_TEMPLATED_INSERTION_LENGTH + 1);
    private static int MAX_TEMPLATED_INSERTION_LENGTH = 500;
    public static int CLUSTER_SIZE_ANALYSIS_LIMIT = 500;
    public static int SHORT_DB_LENGTH = 30;

    public static String TRANS_TYPE_TRANS = "TRANS";
    public static String TRANS_TYPE_SPAN = "SPAN";

    private static final Logger LOGGER = LogManager.getLogger(LinkFinder.class);

    public LinkFinder()
    {
        mLogVerbose = false;
    }

    public void findLinkedPairs(final String sampleId, SvCluster cluster)
    {
        List<SvLinkedPair> assemblyLinkedPairs = createAssemblyLinkedPairs(cluster);
        cluster.setAssemblyLinkedPairs(assemblyLinkedPairs);

        if(cluster.hasLinkingLineElements())
            return;

        List<SvLinkedPair> inferredLinkedPairs = createInferredLinkedPairs(cluster, false);

        // findSpanningSVs(sampleId, cluster, inferredLinkedPairs);

        List<SvLinkedPair> singleBELinkedPairs = createSingleBELinkedPairs(cluster);
        inferredLinkedPairs.addAll(singleBELinkedPairs);

        cluster.setInferredLinkedPairs(inferredLinkedPairs);
    }

    public List<SvLinkedPair> createAssemblyLinkedPairs(SvCluster cluster)
    {
        List<SvLinkedPair> linkedPairs = Lists.newArrayList();

        // find 2 breakends with matching assembly info and form them into a linked pair
        // if have have multiple assembly info and it doesn't match, don't link them
        if(cluster.getCount() < 2)
            return linkedPairs;

        // isSpecificCluster(cluster);

        for (int i = 0; i < cluster.getCount(); ++i)
        {
            SvVarData var1 = cluster.getSVs().get(i);

            if(var1.type() == StructuralVariantType.INS || var1.isNullBreakend())
                continue;

            // make note of SVs which line up exactly with other SVs
            // these will be used to eliminate transitive SVs later on
            if(var1.isDupBreakend(true) && var1.isDupBreakend(false))
            {
                continue;
            }

            for(int be1 = SVI_START; be1 <= SVI_END; ++be1)
            {
                boolean v1Start = isStart(be1);

                for (int j = i+1; j < cluster.getCount(); ++j)
                {
                    SvVarData var2 = cluster.getSVs().get(j);

                    if(var1.equals(var2, true))
                        continue;

                    if(var2.type() == StructuralVariantType.INS || var2.isNullBreakend())
                        continue;

                    if(var2.isDupBreakend(true) && var2.isDupBreakend(false))
                        continue;

                    for(int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        boolean v2Start = isStart(be2);

                        if (!haveLinkedAssemblies(var1, var2, v1Start, v2Start))
                            continue;

                        // check wasn't already created
                        boolean v1Linked = var1.isAssemblyMatched(v1Start);
                        boolean v2Linked = var2.isAssemblyMatched(v2Start);
                        if(v1Linked || v2Linked)
                        {
                            if (v1Linked && v2Linked)
                            {
                                // both linked but to other variants
                            }
                            else if (v1Linked)
                            {
                                var2.setAssemblyMatchType(ASSEMBLY_MATCH_DIFF, v2Start);
                            }
                            else if (v2Linked)
                            {
                                var1.setAssemblyMatchType(ASSEMBLY_MATCH_DIFF, v1Start);
                            }

                            continue;
                        }

                        // form a new TI from these 2 BEs
                        SvLinkedPair newPair = new SvLinkedPair(var1, var2, LINK_TYPE_TI, v1Start, v2Start);
                        newPair.setIsInferred(false);
                        var1.setAssemblyMatchType(ASSEMBLY_MATCH_MATCHED, v1Start);
                        var2.setAssemblyMatchType(ASSEMBLY_MATCH_MATCHED, v2Start);
                        var1.setLinkedPair(newPair, v1Start);
                        var2.setLinkedPair(newPair, v1Start);

                        linkedPairs.add(newPair);

                        // to avoid logging unlikely long TIs
                        LOGGER.debug("cluster({}) adding assembly linked {} pair({}) length({})",
                                cluster.id(), newPair.linkType(), newPair.toString(), newPair.length());
                    }
                }
            }
        }

        return linkedPairs;
    }
    public static boolean areLinkedSection(final SvVarData v1, final SvVarData v2, boolean v1Start, boolean v2Start)
    {
        return areLinkedSection(v1, v2, v1Start, v2Start, false);
    }

    public static boolean areLinkedSection(final SvVarData v1, final SvVarData v2, boolean v1Start, boolean v2Start, boolean checkCopyNumberMatch)
    {
        // templated insertions are allowed to traverse the centromere
        if(v1.position(v1Start) < 0 || v2.position(v1Start) < 0)
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

        if(checkCopyNumberMatch)
        {
            boolean skipReplicated = v1.isReplicatedSv() || v1.getReplicatedCount() > 0 || v2.isReplicatedSv() || v2.getReplicatedCount() > 0;

            if(!skipReplicated)
            {
                double cn1 = v1.copyNumberChange(v1Start);
                double cn2 = v2.copyNumberChange(v2Start);

                if(!copyNumbersEqual(cn1, cn2))
                    return false;
            }
        }

        return true;
    }

    public static boolean areSectionBreak(final SvVarData v1, final SvVarData v2, boolean v1Start, boolean v2Start)
    {
        if(!v1.chromosome(v1Start).equals(v2.chromosome(v2Start)))
            return false;

        // start apart or equal and heading same direction
        long pos1 = v1.position(v1Start);
        boolean headsLeft1 = (v1.orientation(v1Start) == 1);
        long pos2 = v2.position(v2Start);
        boolean headsLeft2 = (v2.orientation(v2Start) == 1);

        if(pos1 <= pos2 && headsLeft1 && !headsLeft2)
            return true;

        if(pos2 <= pos1 && !headsLeft1 && headsLeft2)
            return true;

        return false;
    }


    public List<SvLinkedPair> createInferredLinkedPairs(SvCluster cluster, boolean allowSingleBEs)
    {
        return createInferredLinkedPairs(cluster, cluster.getSVs(), allowSingleBEs);
    }

    public List<SvLinkedPair> createInferredLinkedPairs(SvCluster cluster, List<SvVarData> svList, boolean allowSingleBEs)
    {
        List<SvLinkedPair> linkedPairs = Lists.newArrayList();

        // exclude large clusters for now due to processing times until the algo is better refined
        if(svList.size() >= CLUSTER_SIZE_ANALYSIS_LIMIT)
            return linkedPairs;

        if(cluster.hasLinkingLineElements())
            return linkedPairs;

        for (int i = 0; i < svList.size(); ++i)
        {
            SvVarData var1 = svList.get(i);

            if(var1.type() == StructuralVariantType.INS || (var1.isNullBreakend() && !allowSingleBEs))
                continue;

            if(var1.isDupBreakend(true) && var1.isDupBreakend(false))
                continue;

            for(int be1 = SVI_START; be1 <= SVI_END; ++be1)
            {
                boolean v1Start = isStart(be1);

                if(var1.isNullBreakend() && !v1Start)
                    continue;

                // if an assembly linked pair has already been created for this breakend, look no further
                if(var1.isAssemblyMatched(v1Start))
                    continue;

                for (int j = i+1; j < svList.size(); ++j)
                {
                    SvVarData var2 = svList.get(j);

                    if(var1.equals(var2, true))
                        continue;

                    if(var2.type() == StructuralVariantType.INS || (var2.isNullBreakend() && !allowSingleBEs))
                        continue;

                    if(var2.isDupBreakend(true) && var2.isDupBreakend(false))
                        continue;

                    for(int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        boolean v2Start = isStart(be2);

                        if(var1.isNullBreakend() && !v1Start)
                            continue;

                        if(var2.isAssemblyMatched(v2Start))
                            continue;

                        if (!areLinkedSection(var1, var2, v1Start, v2Start))
                            continue;

                            // form a new TI from these 2 BEs
                        SvLinkedPair newPair = new SvLinkedPair(var1, var2, LINK_TYPE_TI, v1Start, v2Start);

                        if(newPair.linkType() == LINK_TYPE_DB)
                        {
                            // was considered too short to be a TI and so converted
                            var1.setDBLink(newPair, v1Start);
                            var2.setDBLink(newPair, v2Start);
                            continue;
                        }

                        // insert in order of increasing length
                        int index = 0;
                        boolean skipNewPair = false;
                        int linkClashCount = 0;

                        for (; index < linkedPairs.size(); ++index)
                        {
                            final SvLinkedPair pair = linkedPairs.get(index);

                            // check for matching BEs on a pair that is much shorter, and if so skip creating this new linked pair
                            if(newPair.length() > pair.length() && newPair.hasLinkClash(pair))
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
                        if(mLogVerbose && !var1.isReplicatedSv() && !var2.isReplicatedSv())
                        {
                            // to avoid logging unlikely long TIs
                            LOGGER.debug("cluster({}) adding inferred linked {} pair({}) length({}) at index({})",
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

        LOGGER.debug("cluster({}) has {} inferred linked pairs", cluster.id(), linkedPairs.size());

        return linkedPairs;
    }

    public static void findDeletionBridges(final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        for (final Map.Entry<String, List<SvBreakend>> entry : chrBreakendMap.entrySet())
        {
            final List<SvBreakend> breakendList = entry.getValue();

            for (int i = 0; i < breakendList.size() - 1; ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvBreakend nextBreakend = breakendList.get(i+1);
                SvVarData var1 = breakend.getSV();
                SvVarData var2 = nextBreakend.getSV();

                if(var1 == var2)
                    continue;

                isSpecificSV(var1.id());
                isSpecificSV(var2.id());

                if(areSectionBreak(var1, var2, breakend.usesStart(), nextBreakend.usesStart()))
                {
                    SvLinkedPair dbPair = new SvLinkedPair(var1, var2, LINK_TYPE_DB, breakend.usesStart(), nextBreakend.usesStart());
                    markDeletionBridge(dbPair);
                }
                else if(areLinkedSection(var1, var2, breakend.usesStart(), nextBreakend.usesStart())
                && nextBreakend.position() - breakend.position() < MIN_TEMPLATED_INSERTION_LENGTH)
                {
                    // will be converted into a DB
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

    public static boolean arePairedDeletionBridges(final SvVarData var1, final SvVarData var2)
    {
        if(var1.getDBLink(true) == null || var1.getDBLink(false) == null
        || var2.getDBLink(true) == null || var2.getDBLink(false) == null)
        {
            return false;
        }

        if(var1.getDBLink(true) == var2.getDBLink(true) && var1.getDBLink(false) == var2.getDBLink(false))
            return true;
        else if(var1.getDBLink(true) == var2.getDBLink(false) && var1.getDBLink(false) == var2.getDBLink(true))
            return true;
        else
            return false;
    }

    public List<SvLinkedPair> createSingleBELinkedPairs(SvCluster cluster)
    {
        List<SvLinkedPair> linkedPairs = Lists.newArrayList();

        for (int i = 0; i < cluster.getCount(); ++i)
        {
            SvVarData var1 = cluster.getSVs().get(i);

            if(!var1.isNullBreakend() || var1.isNoneSegment())
                continue;

            for (int j = i+1; j < cluster.getCount(); ++j)
            {
                SvVarData var2 = cluster.getSVs().get(j);

                if(var1.equals(var2, true))
                    continue;

                if(!var2.isNullBreakend() || var2.isNoneSegment())
                    continue;

                if (!areSectionBreak(var1, var2, true, true) && !areLinkedSection(var1, var2, true, true))
                {
                    continue;
                }

                // allow if the length is within the DB-TI cutoff
                long length = getProximity(var1, var2, true, true);

                if(length > MIN_TEMPLATED_INSERTION_LENGTH)
                    continue;

                // form a new from these 2 single breakends
                SvLinkedPair newPair = new SvLinkedPair(var1, var2, SvLinkedPair.LINK_TYPE_SGL, false, false);

                // insert in order
                int index = 0;
                boolean skipNewPair = false;
                for (; index < linkedPairs.size(); ++index)
                {
                    SvLinkedPair pair = linkedPairs.get(index);

                    // check for a matching BE on a pair that is shorter, and if so skip creating this new linked pair
                    if(pair.hasAnySameVariant(newPair) && newPair.length() > pair.length())
                    {
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

                if(newPair.length() < DEFAULT_PROXIMITY_DISTANCE)
                {
                    // to avoid logging unlikely long TIs
                    LOGGER.debug("cluster({}) adding inferred single-BE linked {} pair({}) length({}) at index({})",
                            cluster.id(), newPair.linkType(), newPair.toString(), newPair.length(), index);
                }
            }
        }

        if(!linkedPairs.isEmpty())
        {
            LOGGER.debug("cluster({}) has {} inferred single-BE linked pairs",
                    cluster.id(), linkedPairs.size());
        }

        return linkedPairs;
    }

}

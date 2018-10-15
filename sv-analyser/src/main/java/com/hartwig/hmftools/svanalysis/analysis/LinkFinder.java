package com.hartwig.hmftools.svanalysis.analysis;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.haveLinkedAssemblies;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_DIFF;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_DB;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.LINK_TYPE_TI;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class LinkFinder
{
    final SvClusteringConfig mConfig;
    final SvUtilities mUtils;
    SvClusteringMethods mClusteringMethods;

    // private ChainFinder mChainFinder;

    public static int MIN_TEMPLATED_INSERTION_LENGTH = 30;
    private static int MAX_TEMPLATED_INSERTION_LENGTH = 500;
    // public static double MAX_COPY_NUMBER_DIFF = 0.5;
    // public static double MAX_COPY_NUMBER_DIFF_PERC = 0.1;
    public static int CLUSTER_SIZE_ANALYSIS_LIMIT = 200;

    public static String TRANS_TYPE_TRANS = "TRANS";
    public static String TRANS_TYPE_SPAN = "SPAN";

    private static final Logger LOGGER = LogManager.getLogger(LinkFinder.class);

    public LinkFinder(final SvClusteringConfig config, final SvUtilities utils, SvClusteringMethods clusteringMethods)
    {
        mConfig = config;
        mUtils = utils;
        // mChainFinder = new ChainFinder(mUtils);
        // mChainFinder.setLogVerbose(mConfig.LogVerbose);
    }

    public void findLinkedPairs(final String sampleId, SvCluster cluster)
    {
        List<SvLinkedPair> assemblyLinkedPairs = createAssemblyLinkedPairs(sampleId, cluster);
        List<SvLinkedPair> inferredLinkedPairs = createInferredLinkedPairs(sampleId, cluster, false);

        findSpanningSVs(sampleId, cluster, inferredLinkedPairs);

        List<SvLinkedPair> singleBELinkedPairs = createSingleBELinkedPairs(sampleId, cluster);
        inferredLinkedPairs.addAll(singleBELinkedPairs);

        cluster.setAssemblyLinkedPairs(assemblyLinkedPairs);
        cluster.setInferredLinkedPairs(inferredLinkedPairs);
    }

    public List<SvLinkedPair> createAssemblyLinkedPairs(final String sampleId, SvCluster cluster)
    {
        List<SvLinkedPair> linkedPairs = Lists.newArrayList();

        // find 2 breakends with matching assembly info and form them into a linked pair
        // if have have multiple assembly info and it doesn't match, don't link them
        if(cluster.getCount() < 2)
            return linkedPairs;

        for (int i = 0; i < cluster.getCount(); ++i) {

            SvClusterData var1 = cluster.getSVs().get(i);

            if(var1.type() == StructuralVariantType.INS || var1.isNullBreakend())
                continue;

            // make note of SVs which line up exactly with other SVs
            // these will be used to eliminate transitive SVs later on
            if(var1.isDupBEStart() && var1.isDupBEEnd())
            {
                continue;
            }

            for(int be1 = SVI_START; be1 <= SVI_END; ++be1)
            {
                boolean v1Start = isStart(be1);

                for (int j = i+1; j < cluster.getCount(); ++j)
                {
                    SvClusterData var2 = cluster.getSVs().get(j);

                    if(var2.type() == StructuralVariantType.INS || var2.isNullBreakend())
                        continue;

                    if(var2.isDupBEStart() && var2.isDupBEEnd())
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

                        linkedPairs.add(newPair);

                        // to avoid logging unlikely long TIs
                        LOGGER.debug("sample({}) cluster({}) adding assembly linked {} pair({}) length({})",
                                sampleId, cluster.getId(), newPair.linkType(), newPair.toString(), newPair.length());
                    }
                }
            }
        }

        return linkedPairs;
    }

    public List<SvLinkedPair> createInferredLinkedPairs(final String sampleId, SvCluster cluster, boolean allowSingleBEs)
    {
        List<SvLinkedPair> linkedPairs = Lists.newArrayList();

        // exclude large clusters for now due to processing times until the algo is better refined
        if(cluster.getCount() >= CLUSTER_SIZE_ANALYSIS_LIMIT)
            return linkedPairs;

        if(cluster.hasLinkingLineElements())
            return linkedPairs;

        for (int i = 0; i < cluster.getCount(); ++i)
        {
            SvClusterData var1 = cluster.getSVs().get(i);

            if(var1.type() == StructuralVariantType.INS || (var1.isNullBreakend() && !allowSingleBEs))
                continue;

            if(var1.isDupBEStart() && var1.isDupBEEnd())
                continue;

            for(int be1 = SVI_START; be1 <= SVI_END; ++be1)
            {
                boolean v1Start = isStart(be1);

                if(var1.isNullBreakend() && !v1Start)
                    continue;

                // if an assembly linked pair has already been created for this breakend, look no further
                if(var1.isAssemblyMatched(v1Start))
                    continue;

                for (int j = i+1; j < cluster.getCount(); ++j)
                {
                    SvClusterData var2 = cluster.getSVs().get(j);

                    if(var2.type() == StructuralVariantType.INS || (var2.isNullBreakend() && !allowSingleBEs))
                        continue;

                    if(var2.isDupBEStart() && var2.isDupBEEnd())
                        continue;

                    for(int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        boolean v2Start = isStart(be2);

                        if(var1.isNullBreakend() && !v1Start)
                            continue;

                        if(var2.isAssemblyMatched(v2Start))
                            continue;

                        SvLinkedPair newPair = null;

                        if (mUtils.areLinkedSection(var1, var2, v1Start, v2Start))
                        {
                            // form a new TI from these 2 BEs
                            newPair = new SvLinkedPair(var1, var2, LINK_TYPE_TI, v1Start, v2Start);
                        }
                        else if (mUtils.areSectionBreak(var1, var2, v1Start, v2Start))
                        {
                            // form a new DB from these 2 BEs
                            newPair = new SvLinkedPair(var1, var2, LINK_TYPE_DB, v1Start, v2Start);
                        }
                        else
                        {
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
                                if(newPair.linkType() == LINK_TYPE_TI && pair.linkType() == LINK_TYPE_TI)
                                    ++linkClashCount;
                                else if(newPair.linkType() == LINK_TYPE_DB)
                                    ++linkClashCount;

                                if(linkClashCount >= 1)
                                {
                                    skipNewPair = true;
                                    break;
                                }
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

                        // if(linkedPairs.size() > CLUSTER_SIZE_ANALYSIS_LIMIT * 100)
                        //     linkedPairs.remove(linkedPairs.size()-1);

                        // if(newPair.length() < mUtils.getBaseDistance())
                        if(mConfig.LogVerbose && !var1.isReplicatedSv() && !var2.isReplicatedSv())
                        {
                            // to avoid logging unlikely long TIs
                            LOGGER.debug("sample({}) cluster({}) adding inferred linked {} pair({}) length({}) at index({})",
                                    sampleId, cluster.getId(), newPair.linkType(), newPair.toString(), newPair.length(), index);
                        }
                    }
                }
            }
        }

        // prior to consolidating linked pairs, check for duplicate BE in the spanning SVs
        // matchDuplicateBEToLinkedPairs(linkedPairs, spanningSVs);

        if(linkedPairs.isEmpty())
            return linkedPairs;

        LOGGER.debug("sample({}) cluster({}) has {} inferred linked pairs",
                sampleId, cluster.getId(), linkedPairs.size());

        return linkedPairs;
    }

    private void findSpanningSVs(final String sampleId, SvCluster cluster, final List<SvLinkedPair> linkedPairs)
    {
        if(cluster.getCount() >= CLUSTER_SIZE_ANALYSIS_LIMIT)
            return;

        if(cluster.hasLinkingLineElements())
            return;

        List<SvClusterData> spanningSVs = Lists.newArrayList();

        for (int i = 0; i < cluster.getCount(); ++i)
        {
            SvClusterData var1 = cluster.getSVs().get(i);

            if (var1.type() == StructuralVariantType.INS || var1.isNullBreakend())
                continue;

            // make note of SVs which line up exactly with other SVs
            // these will be used to eliminate transitive SVs later on
            if (var1.isDupBEStart() && var1.isDupBEEnd())
            {
                spanningSVs.add(var1);
            }
        }

        if(spanningSVs.isEmpty())
            return;

        // prior to consolidating linked pairs, check for duplicate BE in the spanning SVs
        matchDuplicateBEToLinkedPairs(linkedPairs, spanningSVs);

        LOGGER.debug("sample({}) cluster({}) has {} possible spanning SVs", sampleId, cluster.getId(), spanningSVs.size());

        cluster.setSpanningSVs(spanningSVs);
    }

    public List<SvLinkedPair> createSingleBELinkedPairs(final String sampleId, SvCluster cluster)
    {
        List<SvLinkedPair> linkedPairs = Lists.newArrayList();

        for (int i = 0; i < cluster.getCount(); ++i)
        {
            SvClusterData var1 = cluster.getSVs().get(i);

            if(!var1.isNullBreakend())
                continue;

            for (int j = i+1; j < cluster.getCount(); ++j)
            {
                SvClusterData var2 = cluster.getSVs().get(j);

                if(!var2.isNullBreakend())
                    continue;

                if (!mUtils.areSectionBreak(var1, var2, true, true))
                    continue;

                // form a new DB from these 2 BEs
                SvLinkedPair newPair = new SvLinkedPair(var1, var2, SvLinkedPair.LINK_TYPE_SGL, false, false);

                // insert in order
                int index = 0;
                boolean skipNewPair = false;
                for (; index < linkedPairs.size(); ++index)
                {
                    SvLinkedPair pair = linkedPairs.get(index);

                    // check for a matching BE on a pair that is much shorter, and if so skip creating this new linked pair
                    if(newPair.length() > mUtils.getBaseDistance())
                    {
                        if (pair.first().equals(newPair.first()) || pair.first().equals(newPair.second()) || pair.second().equals(newPair.first()) || pair.second().equals(newPair.second())) {

                            if (newPair.length() > 2 * pair.length())
                            {
                                skipNewPair = true;
                                break;
                            }
                        }
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

                if(newPair.length() < mUtils.getBaseDistance())
                {
                    // to avoid logging unlikely long TIs
                    LOGGER.debug("sample({}) cluster({}) adding inferred single-BE linked {} pair({}) length({}) at index({})",
                            sampleId, cluster.getId(), newPair.linkType(), newPair.toString(), newPair.length(), index);
                }
            }
        }

        if(!linkedPairs.isEmpty())
        {
            LOGGER.debug("sample({}) cluster({}) has {} inferred single-BE linked pairs",
                    sampleId, cluster.getId(), linkedPairs.size());
        }

        return linkedPairs;
    }

    private void matchDuplicateBEToLinkedPairs(final List<SvLinkedPair> linkedPairs, final List<SvClusterData> spanningSVs)
    {
        // link spanning SVs with any single linked pairs
        for(SvClusterData spanningSV : spanningSVs)
        {
            SvClusterData startLink = null;
            SvClusterData endLink = null;

            for(SvLinkedPair pair : linkedPairs) {

                if (pair.length() > MAX_TEMPLATED_INSERTION_LENGTH) {
                    continue;
                }

                if(mUtils.breakendsMatch(spanningSV, pair.first(), true, !pair.firstLinkOnStart(), PERMITED_DUP_BE_DISTANCE)) {
                    startLink = pair.first();
                }
                else if(mUtils.breakendsMatch(spanningSV, pair.second(), true, !pair.secondLinkOnStart(), PERMITED_DUP_BE_DISTANCE)) {
                    startLink = pair.second();
                }
                else
                {
                    continue;
                }

                if(mUtils.breakendsMatch(spanningSV, pair.first(), false, !pair.firstLinkOnStart(), PERMITED_DUP_BE_DISTANCE)) {
                    endLink = pair.first();
                }
                else if(mUtils.breakendsMatch(spanningSV, pair.second(), false, !pair.secondLinkOnStart(), PERMITED_DUP_BE_DISTANCE))
                {
                    endLink = pair.second();
                }
                else
                {
                    continue;
                }

                // match found on both ends
                LOGGER.debug("spanSV({}) linked to linked pair({} and {})",
                        spanningSV.posId(), startLink.posId(), endLink.posId());

                startLink.setTransData(TRANS_TYPE_TRANS, pair.length(), spanningSV.id());
                endLink.setTransData(TRANS_TYPE_TRANS, pair.length(), spanningSV.id());

                String svLinkData = startLink.id() + "_" + endLink.id();
                spanningSV.setTransData(TRANS_TYPE_SPAN, pair.length(), svLinkData);

                break;
            }
        }
    }

    public void resolveTransitiveSVs(final String sampleId, SvCluster cluster)
    {
        if (cluster.getLinkedPairs().isEmpty() || cluster.getSpanningSVs().isEmpty())
            return;

        // attempt to matching spanning SVs to the ends of one or more linked pairs
        // these can only span short TIs (ie not DBs or long TIs)
        for(SvClusterData spanningSV : cluster.getSpanningSVs())
        {
            boolean startMatched = false;
            boolean endMatched = false;
            SvClusterData startLink = null;
            SvClusterData endLink = null;
            SvLinkedPair startPair = null;
            SvLinkedPair endPair = null;
            boolean startLinkOnStart = false;
            boolean endLinkOnStart = false;

            List<SvClusterData> transitiveSVs = Lists.newArrayList();

            for(SvLinkedPair pair : cluster.getLinkedPairs()) {

                if (pair.length() > MAX_TEMPLATED_INSERTION_LENGTH) {
                    continue;
                }

                if (!startMatched)
                {
                    if(mUtils.breakendsMatch(spanningSV, pair.first(), true, !pair.firstLinkOnStart(), PERMITED_DUP_BE_DISTANCE)) {
                        startLink = pair.first();
                        startLinkOnStart = !pair.firstLinkOnStart();
                    }
                    else if(mUtils.breakendsMatch(spanningSV, pair.second(), true, !pair.secondLinkOnStart(), PERMITED_DUP_BE_DISTANCE)) {
                        startLink = pair.second();
                        startLinkOnStart = !pair.secondLinkOnStart();
                    }
                    else
                    {
                        continue;
                    }

                    startMatched = true;
                    startPair = pair;
                }

                if (!endMatched)
                {
                    if(mUtils.breakendsMatch(spanningSV, pair.first(), false, !pair.firstLinkOnStart(), PERMITED_DUP_BE_DISTANCE)) {
                        endLink = pair.first();
                        endLinkOnStart = !pair.firstLinkOnStart();
                    }
                    else if(mUtils.breakendsMatch(spanningSV, pair.second(), false, !pair.secondLinkOnStart(), PERMITED_DUP_BE_DISTANCE))
                    {
                        endLink = pair.second();
                        endLinkOnStart = !pair.secondLinkOnStart();
                    }
                    else
                    {
                        continue;
                    }

                    endPair = pair;
                    endMatched = true;
                }


                if(startMatched && endMatched)
                    break;
            }

            if(startMatched && endMatched)
            {
                boolean samePair = (startPair == endPair);
                int tiLength = 0;
                boolean hasValidTransData = true;

                LOGGER.debug("cluster({}) spanSV({}) linked to transitives({} and {}) from {} linked pair",
                        cluster.getId(), spanningSV.posId(), startLink.posId(), endLink.posId(),
                        samePair ? "same" : "diff");

                if(samePair)
                {
                    tiLength = startPair.length();
                    transitiveSVs.add(startLink);
                    transitiveSVs.add(endLink);
                }
                else {

                    // now additionally check if these SVs are part of a chain, and if so whether any intermediary transitive SVs
                    // are also covered by this span
                    SvChain startChain = cluster.findChain(startPair);
                    SvChain endChain = cluster.findChain(endPair);

                    if (startChain != null && endChain == startChain) {

                        // now walk the chain and collect up all transitive SVs
                        int startIndex = startChain.getSvIndex(startLink, startLinkOnStart);
                        int endIndex = startChain.getSvIndex(endLink, endLinkOnStart);

                        if(startIndex == -1 || endIndex == -1)
                        {
                            LOGGER.error("cluster({}) chain({}) index not found({} - {}) links({} & {})",
                                    cluster.getId(), startChain.getId(), startIndex, endIndex, startLink, endLink);
                            return;
                        }

                        int totalTILength = 0;

                        if (startIndex != endIndex) {

                            int startI = startIndex <= endIndex ? startIndex : endIndex;
                            int endI = startIndex > endIndex ? startIndex : endIndex;

                            for(int i = startI; i <= endI; ++i)
                            {
                                final SvLinkedPair pair = startChain.getLinkedPairs().get(i);

                                if(transitiveSVs.contains(pair.first()) || transitiveSVs.contains(pair.second()))
                                {
                                    LOGGER.debug("cluster({}) chain({}) attempt to re-add trans SVs, invalid",
                                            cluster.getId(), startChain.getId());

                                    transitiveSVs.clear();

                                    // manually add the link SVs
                                    totalTILength = 0;
                                    transitiveSVs.add(startLink);
                                    transitiveSVs.add(endLink);
                                    break;
                                }

                                totalTILength += pair.length();

                                if(totalTILength > MAX_TEMPLATED_INSERTION_LENGTH * 3)
                                {
                                    LOGGER.debug("cluster({}) chain({}) exceed valid totalLen({}) at index({}), invalid",
                                            cluster.getId(), startChain.getId(), totalTILength, i);

                                    hasValidTransData = false;
                                    break;
                                }

                                LOGGER.debug("cluster({}) chain({}) including index({}) totalLen({}) linkPair({}))",
                                        cluster.getId(), startChain.getId(), i, totalTILength, pair.toString());

                                transitiveSVs.add(pair.first());
                                transitiveSVs.add(pair.second());
                            }

                            if(hasValidTransData)
                            {
                                LOGGER.info("cluster({}) spanSV({}) covers {} linked pairs",
                                        cluster.getId(), spanningSV.id(), transitiveSVs.size()/2);

                                tiLength = totalTILength;
                            }
                        }
                        else
                        {
                            LOGGER.warn("cluster({}) chain({}) linked pairs have same index({}) but diff linked pair",
                                    cluster.getId(), startChain.getId(), startIndex);
                        }
                    }
                    else if (startChain == null || endChain == null) {

                        // ignore any intermediary linked SVs from the single chain for now
                        tiLength = startPair.length() + endPair.length();

                        if(tiLength < MAX_TEMPLATED_INSERTION_LENGTH * 3) {

                            transitiveSVs.add(startLink);
                            transitiveSVs.add(endLink);
                        }
                        else
                        {
                            hasValidTransData = false;
                        }
                    }
                    else
                    {
                        hasValidTransData = false;
                        LOGGER.debug("cluster({}) linked pairs have diff chains({} and {})",
                                cluster.getId(), startChain.getId(), endChain.getId());
                    }
                }

                if(hasValidTransData) {

                    // mark all transitive SVs
                    for (SvClusterData transSv : transitiveSVs) {
                        String svLinkData = spanningSV.id();
                        transSv.setTransData(TRANS_TYPE_TRANS, tiLength, svLinkData);
                    }

                    // and mark the spanning SV
                    String svLinkData = startLink.id() + "_" + endLink.id();
                    spanningSV.setTransData(TRANS_TYPE_SPAN, tiLength, svLinkData);
                }
            }
        }
    }

}

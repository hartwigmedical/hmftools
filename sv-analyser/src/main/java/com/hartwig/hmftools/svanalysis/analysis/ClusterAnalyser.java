package com.hartwig.hmftools.svanalysis.analysis;

import static java.lang.Math.round;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.PERMITED_DUP_BE_DISTANCE;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_END;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.SVI_START;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.haveLinkedAssemblies;
import static com.hartwig.hmftools.svanalysis.types.SvClusterData.isStart;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_DIFF;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_INFER_ONLY;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_MATCHED;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvCNData;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ClusterAnalyser {

    final SvClusteringConfig mConfig;
    final SvUtilities mUtils;
    SvClusteringMethods mClusteringMethods;

    private ChainFinder mChainFinder;

    public static int MIN_TEMPLATED_INSERTION_LENGTH = 30;
    private static int MAX_TEMPLATED_INSERTION_LENGTH = 500;
    public static double MAX_COPY_NUMBER_DIFF = 0.5;
    public static double MAX_COPY_NUMBER_DIFF_PERC = 0.1;
    public static int CLUSTER_SIZE_ANALYSIS_LIMIT = 100;

    public static String TRANS_TYPE_TRANS = "TRANS";
    public static String TRANS_TYPE_SPAN = "SPAN";

    private static final Logger LOGGER = LogManager.getLogger(ClusterAnalyser.class);

    public ClusterAnalyser(final SvClusteringConfig config, final SvUtilities utils, SvClusteringMethods clusteringMethods)
    {
        mConfig = config;
        mUtils = utils;
        mClusteringMethods = clusteringMethods;
        mChainFinder = new ChainFinder(mUtils);
    }

    public void setClusterStats(SvCluster cluster)
    {
        cluster.setConsistencyCount();

        // record an expression of the types of SVs in this cluster
        cluster.setDesc(cluster.getClusterTypesAsString());

        if(cluster.getCount() > 1)
        {
            LOGGER.debug("cluster({}) svCount({}) desc({}) armCount({}) consistent({} count={})",
                    cluster.getId(), cluster.getCount(), cluster.getDesc(),
                    cluster.getChromosomalArmCount(), cluster.isConsistent(), cluster.getConsistencyCount());
        }
    }

    public void findLinksAndChains(final String sampleId, List<SvCluster> clusters)
    {
        for(SvCluster cluster : clusters)
        {
            applyCopyNumberReplication(sampleId, cluster);

            if(cluster.getCount() < 2)
                continue;

            // first establish links between SVs (eg TIs and DBs)
            findLinkedPairs(sampleId, cluster);

            // then look for fully-linked clusters, ie chains involving all SVs
            findCompleteChains(sampleId, cluster);
        }

        // now look at merging unresolved & inconsistent clusters where they share the same chromosomal arms
        mergeInconsistentClusters(clusters);

        for(SvCluster cluster : clusters)
        {
            if(cluster.hasSubClusters())
            {
                cluster.setDesc(cluster.getClusterTypesAsString());

                // repeat the search for inferred links now that additional SVs have been merged in
                List<SvLinkedPair> newLinkedPairs = createInferredLinkedPairs(sampleId, cluster, true);
                cluster.getInferredLinkedPairs().addAll(newLinkedPairs);

                createCopyNumberSegments(sampleId, cluster);

                findIncompleteChains(sampleId, cluster);
            }

            resolveTransitiveSVs(sampleId, cluster);
            cacheFinalLinkedPairs(sampleId, cluster);
        }
    }

    private void findLinkedPairs(final String sampleId, SvCluster cluster)
    {
        List<SvLinkedPair> assemblyLinkedPairs = createAssemblyLinkedPairs(sampleId, cluster);
        List<SvLinkedPair> inferredLinkedPairs = createInferredLinkedPairs(sampleId, cluster, false);

        findSpanningSVs(sampleId, cluster, inferredLinkedPairs);

        List<SvLinkedPair> singleBELinkedPairs = createSingleBELinkedPairs(sampleId, cluster);
        inferredLinkedPairs.addAll(singleBELinkedPairs);

        cluster.setAssemblyLinkedPairs(assemblyLinkedPairs);
        cluster.setInferredLinkedPairs(inferredLinkedPairs);
        cluster.setConsistencyCount(); // since assembly replicated links can change consistency
    }

    private void findCompleteChains(final String sampleId, SvCluster cluster)
    {
        mChainFinder.initialise(sampleId, cluster);
        mChainFinder.setRequireFullChains(true);

        boolean hasFullChain = mChainFinder.formClusterChains();

        if(!hasFullChain)
            return;

        cluster.setIsFullyChained(true);

        // remove any inferred link which isn't in the full chain
        final SvChain fullChain = cluster.getChains().get(0);

        List<SvLinkedPair> inferredLinkedPairs = cluster.getInferredLinkedPairs();

        // remove any inferred links which weren't used
        int lpIndex = 0;
        while(lpIndex < inferredLinkedPairs.size())
        {
            SvLinkedPair pair = inferredLinkedPairs.get(lpIndex);

            if (!fullChain.getLinkedPairs().contains(pair))
            {
                inferredLinkedPairs.remove(lpIndex);
                continue;
            }

            ++lpIndex;
        }
    }

    private void applyCopyNumberReplication(final String sampleId, SvCluster cluster)
    {
        List<SvClusterData> newSVs = Lists.newArrayList();

        for(final SvClusterData var : cluster.getSVs())
        {
            int calcCopyNumber = var.impliedCopyNumber(true);

            int chromosomeCopyNumber = mClusteringMethods.getChrCopyNumberMap().get(var.chromosome(true));

            // TEMP: until chromosomal copy-number calc working
            if(mConfig.SampleCopyNumber != 1)
            {
                chromosomeCopyNumber = mConfig.SampleCopyNumber;
            }

            if(chromosomeCopyNumber >= 2)
                calcCopyNumber /= chromosomeCopyNumber;

            if(calcCopyNumber <= 1)
                continue;

            LOGGER.debug("sample({}) replicating SV({}) {} times", sampleId, var.posId(), calcCopyNumber-1);

            for(int i = 1; i < calcCopyNumber; ++i)
            {
                SvClusterData newVar = new SvClusterData(var);
                newSVs.add(newVar);
            }
        }

        for(SvClusterData var : newSVs)
        {
            cluster.addVariant(var);
        }
    }

    private void mergeInconsistentClusters(List<SvCluster> clusters)
    {
        // it's possible that to resolve arms and more complex arrangements, that clusters not merged
        // by proximity of overlaps must be put together to solve inconsistencies (ie loose ends)

        // merge any cluster which is itself not consistent and has breakends on the same arm as another
        int index1 = 0;
        while(index1 < clusters.size())
        {
            SvCluster cluster1 = clusters.get(index1);

            if(cluster1.isConsistent())
            {
                ++index1;
                continue;
            }

            if(!cluster1.isFullyChained() && cluster1.getUniqueSvCount() > 1 && !cluster1.hasSubClusters())
            {
                ++index1;
                continue;
            }

            boolean cluster1Merged = false;
            SvCluster newCluster = null;

            int index2 = index1 + 1;
            while(index2 < clusters.size())
            {
                SvCluster cluster2 = clusters.get(index2);

                if(cluster2.isConsistent())
                {
                    ++index2;
                    continue;
                }

                if(!cluster2.isFullyChained() && cluster2.getUniqueSvCount() > 1 && !cluster1.hasSubClusters())
                {
                    ++index2;
                    continue;
                }

                boolean cluster2Merged = false;

                boolean foundConnection = false;

                for(int be1 = SVI_START; be1 <= SVI_END; ++be1)
                {
                    final String c1LinkingArm = cluster1.linkingArm(isStart(be1));

                    if(c1LinkingArm.isEmpty())
                        continue;

                    final String c1LinkingChr = cluster1.linkingChromosome(isStart(be1));

                    for(int be2 = SVI_START; be2 <= SVI_END; ++be2)
                    {
                        final String c2LinkingArm = cluster2.linkingArm(isStart(be2));

                        if(c2LinkingArm.isEmpty())
                            continue;

                        if(c1LinkingArm.equals(c2LinkingArm) && c1LinkingChr.equals(cluster2.linkingChromosome(isStart(be2))))
                        {
                            LOGGER.debug("cluster({}) and cluster({}) linked on chrArm({}:{})",
                                    cluster1.getId(), cluster2.getId(), c1LinkingChr, c1LinkingArm);

                            foundConnection = true;
                            break;
                        }
                    }

                    if(foundConnection)
                        break;
                }

                if(!foundConnection)
                {
                    ++index2;
                    continue;
                }

                if(cluster1.hasSubClusters())
                {
                    LOGGER.debug("cluster({} svs={}) merges in cluster({} svs={})",
                            cluster1.getId(), cluster1.getCount(), cluster2.getId(), cluster2.getCount());

                    cluster2Merged = true;
                    cluster1.addSubCluster(cluster2);
                }
                else
                {
                    cluster1Merged = true;
                    cluster2Merged = true;

                    newCluster = new SvCluster(mClusteringMethods.getNextClusterId(), mUtils);

                    LOGGER.debug("new cluster({}) from merge of cluster({} svs={}) and cluster({} svs={})",
                            newCluster.getId(), cluster1.getId(), cluster1.getCount(), cluster2.getId(), cluster2.getCount());

                    newCluster.addSubCluster(cluster1);
                    newCluster.addSubCluster(cluster2);
                }

                if(cluster2Merged)
                    clusters.remove(index2);
                else
                    ++index2;

                if(cluster1Merged)
                    break;
            }

            if(cluster1Merged && newCluster != null)
            {
                // cluster has been replaced with a commbined one
                clusters.remove(index1);
                clusters.add(index1, newCluster);
            }
            else
            {
                ++index1;
            }
        }
    }

    private void findIncompleteChains(final String sampleId, SvCluster cluster)
    {
        if(cluster.getCount() < 2 || cluster.isFullyChained())
            return;

        mChainFinder.initialise(sampleId, cluster);
        mChainFinder.setRequireFullChains(false);

        mChainFinder.formClusterChains();
    }

    private void createCopyNumberSegments(final String sampleId, final SvCluster cluster)
    {
        int cnId = 0;
        final List<SvClusterData> clusterSVs = cluster.getSVs();

        Map<String, List<SvCNData>> chrCNDataMap = new HashMap();

        for(Map.Entry<String, List<SvBreakend>> entry : mClusteringMethods.getChrBreakendMap().entrySet())
        {
            final String chromosome = entry.getKey();

            List<SvBreakend> breakendList = entry.getValue();

            List<SvCNData> copyNumberData = Lists.newArrayList();

            for (int i = 0; i < breakendList.size(); ++i)
            {
                final SvBreakend breakend = breakendList.get(i);
                final SvClusterData var = breakend.getSV();

                double copyNumber = var.copyNumber(breakend.usesStart());
                double copyNumberChange = var.copyNumberChange(breakend.usesStart());
                double prevCopyNumber = copyNumber - copyNumberChange;

                int adjCopyNumber = (int)round(copyNumber/2 - 1);

                //LOGGER.debug(String.format("sample(%s) chr(%s) seg %d: copyNumber(%d %.2f prev=%.2f chg=%.2f)",
                //        sampleId, chromosome, i, adjCopyNumber, copyNumber, prevCopyNumber, copyNumberChange));

                if(clusterSVs.contains(var))
                {
                    SvCNData cnData = new SvCNData(cnId++, chromosome, breakend.position(), 0, "", "", 0, 0, 0, adjCopyNumber, "");
                    copyNumberData.add(cnData);
                }
            }

            chrCNDataMap.put(chromosome, copyNumberData);
        }

        cluster.setChrCNData(chrCNDataMap);
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
                        boolean allowDuplicateLink = false;

                        /*

                        handled by replicated SVs instead..

                        if(v1Linked || v2Linked)
                        {
                            // check for special case of A-B and A-C where B and C are the same variant (typically a fold-back replicated inversion)
                            // and in this case duplicate the breakend of A so a new link is created
                            if (!v1Linked && v2Linked)
                            {
                                // check if there is already a link from v1's unlinked end to v2's linked end
                                if (haveLinkedAssemblies(var1, var2, !v1Start, v2Start))
                                {
                                    var2.setIsReplicatedLink(!v2Start, true);
                                    allowDuplicateLink = true;
                                }
                            }
                            else if (v1Linked && !v2Linked)
                            {
                                // check if there is already a link from v2's unlinked end to v1's linked end
                                if (haveLinkedAssemblies(var1, var2, v1Start, !v2Start))
                                {
                                    var1.setIsReplicatedLink(!v1Start, true);
                                    allowDuplicateLink = true;
                                }
                            }
                        }
                        */

                        if((v1Linked || v2Linked) && !allowDuplicateLink)
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
                        SvLinkedPair newPair = new SvLinkedPair(var1, var2, SvLinkedPair.LINK_TYPE_TI, v1Start, v2Start);
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

        if(cluster.getLineElementCount() > 0)
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
                            newPair = new SvLinkedPair(var1, var2, SvLinkedPair.LINK_TYPE_TI, v1Start, v2Start);
                        }
                        else if (mUtils.areSectionBreak(var1, var2, v1Start, v2Start))
                        {
                            // form a new DB from these 2 BEs
                            newPair = new SvLinkedPair(var1, var2, SvLinkedPair.LINK_TYPE_DB, v1Start, v2Start);
                        }
                        else
                        {
                            continue;
                        }

                        // insert in order
                        int index = 0;
                        boolean skipNewPair = false;
                        for (; index < linkedPairs.size(); ++index)
                        {
                            SvLinkedPair pair = linkedPairs.get(index);

                            /*
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
                            */

                            if (pair.length() > newPair.length())
                                break;
                        }

                        if(skipNewPair)
                            continue;

                        if (index >= linkedPairs.size())
                            linkedPairs.add(newPair);
                        else
                            linkedPairs.add(index, newPair);

                        if(linkedPairs.size() > CLUSTER_SIZE_ANALYSIS_LIMIT * 3)
                            linkedPairs.remove(linkedPairs.size()-1);

                        if(newPair.length() < mUtils.getBaseDistance())
                        {
                            // to avoid logging unlikely long TIs
                            //LOGGER.debug("sample({}) cluster({}) adding inferred linked {} pair({}) length({}) at index({})",
                            //         sampleId, cluster.getId(), newPair.linkType(), newPair.toString(), newPair.length(), index);
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

        if(cluster.getLineElementCount() > 0)
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

    private void cacheFinalLinkedPairs(final String sampleId, SvCluster cluster)
    {
        List<SvLinkedPair> linkedPairs;

        if(cluster.isFullyChained())
        {
            linkedPairs = cluster.getChains().get(0).getLinkedPairs();
        }
        else
        {
            List<SvLinkedPair> inferredLinkedPairs = cluster.getInferredLinkedPairs();

            reduceInferredToShortestLinks(inferredLinkedPairs, cluster.getChains());

            List<SvLinkedPair> assemblyLinkedPairs = cluster.getAssemblyLinkedPairs();

            // remove any linked which aren't part of chains

            // mark the resultant set of inferred links
            for (SvLinkedPair pair : inferredLinkedPairs)
            {
                pair.first().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.firstLinkOnStart());
                pair.second().setAssemblyMatchType(ASSEMBLY_MATCH_INFER_ONLY, pair.secondLinkOnStart());
            }

            linkedPairs = Lists.newArrayList();
            linkedPairs.addAll(assemblyLinkedPairs);
            linkedPairs.addAll(inferredLinkedPairs);
        }

        if(!linkedPairs.isEmpty())
        {
            LOGGER.info("sample({}) cluster({}: {} count={}) has {} linked pairs:",
                    sampleId, cluster.getId(), cluster.getDesc(), cluster.getUniqueSvCount(), linkedPairs.size());

            for (final SvLinkedPair pair : linkedPairs)
            {
                LOGGER.info("linked {} {} pair length({}) variants({})",
                        pair.isInferred() ? "inferred" : "assembly", pair.linkType(), pair.length(), pair.toString());
            }
        }

        cluster.setLinkedPairs(linkedPairs);
    }

    public static void reduceInferredToShortestLinks(List<SvLinkedPair> linkedPairs, List<SvChain> chains)
    {
        // any linked pair used in a chain must be kept
        List<SvLinkedPair> reqLinkedPairs = Lists.newArrayList();
        for(SvChain chain : chains)
        {
            reqLinkedPairs.addAll(chain.getLinkedPairs());
        }

        // now remove mutually exclusive linked sections by using the shortest first (they are already ordered)
        int i = 0;
        while(i < linkedPairs.size())
        {
            final SvLinkedPair pair = linkedPairs.get(i);

            if(reqLinkedPairs.contains(pair))
            {
                ++i;
                continue;
            }

            boolean removeFirst = false;

            // search for another pair with a matching breakend
            for(int j = i+1; j < linkedPairs.size();)
            {
                final SvLinkedPair pair2 = linkedPairs.get(j);

                if((pair.first().equals(pair2.first()) && pair.firstLinkOnStart() == pair2.firstLinkOnStart())
                || (pair.first().equals(pair2.second()) && pair.firstLinkOnStart() == pair2.secondLinkOnStart())
                || (pair.second().equals(pair2.first()) && pair.secondLinkOnStart() == pair2.firstLinkOnStart())
                || (pair.second().equals(pair2.second()) && pair.secondLinkOnStart() == pair2.secondLinkOnStart()))
                {
                    // to avoid logging unlikely long TIs
                    //LOGGER.debug("removing duplicate linked pair({} len={}) vs shorter({} len={})",
                    //        pair2.toString(), pair2.length(), pair.toString(), pair.length());

                    if(reqLinkedPairs.contains(pair2))
                    {
                        removeFirst = true;
                        break;
                    }

                    // remove the longer pair with a duplicate breakend
                    linkedPairs.remove(j);
                }
                else
                {
                    ++j;
                }
            }

            if(removeFirst)
            {
                // LOGGER.debug("duplicate linked pair({} len={})", pair.toString(), pair.length(), pair.first().getTransSvLinks());
                linkedPairs.remove(i);
            }
            else
            {
                ++i;
            }
        }
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
                        LOGGER.info("cluster({}) linked pairs have diff chains({} and {})",
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

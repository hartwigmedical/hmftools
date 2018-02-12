package com.hartwig.hmftools.svannotation.analysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.svannotation.analysis.SvCluster;
import com.hartwig.hmftools.svannotation.analysis.SvClusteringConfig;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class StructuralVariantClustering {

    private final SvClusteringConfig mConfig;
    private final SvUtilities mClusteringUtils;

    // data per run (ie sample)
    private String mSampleId;
    private List<SvClusterData> mAllVariants; // the original list to analyse
    private List<SvClusterData> mUnassignedVariants; // the working list to be assigned
    private int mNextClusterId;

    private List<SvCluster> mClusters;

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnalyzer.class);

    public StructuralVariantClustering(final SvClusteringConfig config)
    {
        mConfig = config;
        mClusteringUtils = new SvUtilities(mConfig.getClusterBaseDistance());
        clearState();
    }

    private void clearState()
    {
        mSampleId = "";
        mNextClusterId = 0;
        mAllVariants = Lists.newArrayList();
        mUnassignedVariants = Lists.newArrayList();
        mClusters = Lists.newArrayList();
    }

    public void loadFromEnrichedSVs(final String sampleId, final List<EnrichedStructuralVariant> variants)
    {
        if (variants.isEmpty())
            return;

        clearState();

        mSampleId = sampleId;

        for (final EnrichedStructuralVariant enrichedSV : variants)
        {
            mUnassignedVariants.add(SvClusterData.from(enrichedSV));
        }

        mAllVariants = Lists.newArrayList(mUnassignedVariants);
    }

    public void loadFromDatabase(final String sampleId, final List<SvClusterData> variants)
    {
        if (variants.isEmpty())
            return;

        clearState();

        mSampleId = sampleId;
        mUnassignedVariants = Lists.newArrayList(variants);
        mAllVariants = Lists.newArrayList(mUnassignedVariants);

        LOGGER.debug("loaded {} SVs", mAllVariants.size());

//        for(SvClusterData sv : variants) {
//            LOGGER.debug("adding SV: id({}) start({}:{}) end({}:{})",
//                    sv.id(), sv.chromosome(true), sv.position(true), sv.chromosome(false), sv.position(false));
//        }
    }

    public void runClustering()
    {
        if(mUnassignedVariants.isEmpty())
            return;

        LOGGER.debug("sample({}) clustering {} variants", mSampleId, mUnassignedVariants.size());
        // original runClusteringByBaseDistance

        runClusteringByLocals();
    }

    public void runClusteringByLocals() {
        // first find outer-most clusters at CRMS level and identify others which will span CRMS

        // LOGGER.debug("sample({}) isolating inter-chromosomal variants", mSampleId);

        // first remove cross-chromosomal from initial consideration
        // List<SvClusterData> crossChromosomal = findCrossChromosomalSVs();

        LOGGER.debug("sample({}) clustering contained variants", mSampleId);

        createOutermostChromosomalClusters();

        LOGGER.debug("sample({}) assigning sub-clusters", mSampleId);

        // now need to assign and create sub-clusters within these outer-most clusters
        for(SvCluster cluster : mClusters) {
            assignSubClusters(cluster);
        }

        LOGGER.debug("sample({}) matching {} overlapping variants", mSampleId, mUnassignedVariants.size());

        // final create links for the cross-chromosomal SVs and other overlapping SVs
        assignCrossClusterSVs();
    }

    private List<SvClusterData> findCrossChromosomalSVs()
    {
        // first remove cross-chromosomal from initial consideration
        List<SvClusterData> crossChromosomal = Lists.newArrayList();

        int currentIndex = 0;
        while (currentIndex < mUnassignedVariants.size()) {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            // first skip over cross-CRMS for now
            if (!currentVar.isLocal()) {
                crossChromosomal.add(currentVar);

                LOGGER.debug("cross-CRMS SV: {}", currentVar.posId());

                mUnassignedVariants.remove(currentIndex); // index will point at next
                continue;
            } else {
                ++currentIndex;
            }
        }

        return crossChromosomal;

    }

    private void createOutermostChromosomalClusters()
    {
        int currentIndex = 0;

        List<SvClusterData> subSVsToRemove = Lists.newArrayList();

        while(currentIndex < mUnassignedVariants.size()) {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            if(!currentVar.isLocal())
            {
                LOGGER.debug("skipping inter-chromosomal SV: {}", currentVar.posId());
                ++currentIndex;
                continue;
            }

            boolean spansOtherSVs = false;
            List<SvClusterData> subSVs = Lists.newArrayList();

            for (SvClusterData otherVar : mUnassignedVariants) 
            {
                if(otherVar.id().equals(currentVar.id()))
                {
                    continue;
                }

                if(mClusteringUtils.isLocalOverlap(currentVar, otherVar))
                {
                    LOGGER.debug("local overlap SVs: v1({}) and v2({})", currentVar.posId(), otherVar.posId());
                    spansOtherSVs = true;
                }

                if (mClusteringUtils.isWithin(currentVar, otherVar))
                {
                    if(currentVar.addSubSV(otherVar)) {
                        LOGGER.debug("wholy contained: outer({}) and inner({})", currentVar.posId(), otherVar.posId());
                    }

                }
                else if (mClusteringUtils.isWithin(otherVar, currentVar))
                {
                    if(otherVar.addSubSV(currentVar)) {
                        LOGGER.debug("wholy contained: outer({}) and inner({})", otherVar.posId(), currentVar.posId());
                    }
                }
            }

            if(!currentVar.isSubSV() && (currentVar.hasSubSVs() || !spansOtherSVs))
            {
                LOGGER.debug("adding outer-most CRMS SV: {}", currentVar.posId());

                // make a new cluster
                SvCluster newCluster = new SvCluster(getNextClusterId(), mClusteringUtils);

                newCluster.setSpanningSV(currentVar);
                mClusters.add(newCluster);

                mUnassignedVariants.remove(currentIndex); // index will point at next
            }
            else
            {
                if(currentVar.isSubSV())
                    subSVsToRemove.add(currentVar);

                ++currentIndex;
            }
        }

        for(SvClusterData variant : subSVsToRemove)
        {
            // remove from unassigned
            mUnassignedVariants.remove(variant);
        }
    }

    private void assignSubClusters(SvCluster cluster)
    {
        LOGGER.debug("cluster({}) sub-clustering {} variants", cluster.getClusterId(), cluster.getSVs().size());

        // create sub-clusters for any SVs which have sub SVs
        for(SvClusterData currentVar : cluster.getSVs())
        {
            if(currentVar == cluster.getSpanningSV())
                continue;

            if(!currentVar.hasSubSVs())
            {
                if(currentVar.areClustersSet())
                {
                    // keep with the smallest cluster
                    long clusterSpan = currentVar.getStartCluster().getSpan();
                    long newSpan = cluster.getSpan();

                    if(newSpan == -1 || newSpan > clusterSpan)
                    {
                        continue;
                    }

                }

                LOGGER.debug("variant({}) assigned to cluster({}) span({} -> {})",
                        currentVar.posId(), cluster.getClusterId(), cluster.getSpanningSV().position(true), cluster.getSpanningSV().position(false));

                currentVar.setStartCluster(cluster);
                currentVar.setEndCluster(cluster);
                continue;
            }
            else {
                // also check that this variant is a sub-variant of another variant at this level
                boolean skipClustering = false;
                for(SvClusterData other : cluster.getSVs())
                {
                    if(currentVar.equals(other) || other == cluster.getSpanningSV())
                        continue;

                    if(other.getSubSVs().contains(currentVar))
                    {
                        LOGGER.debug("skipping sub-cluster for variant({}) since part of other({})", currentVar.id(), other.id());
                        skipClustering = true;
                        break;
                    }
                }

                if(skipClustering)
                    continue;
            }

            // create a cluster, set its parent cluster and then call iteratively
            SvCluster newCluster = new SvCluster(getNextClusterId(), mClusteringUtils);

            LOGGER.debug("cluster({}) creating new subCluster({})", cluster.getClusterId(), newCluster.getClusterId());

            // add the spanning SV and its sub-SVs
            newCluster.setSpanningSV(currentVar);
            cluster.addSubCluster(newCluster);

            // and call recursively to create lower-level clusters
            assignSubClusters(newCluster);
        }
    }

    private void assignCrossClusterSVs()
    {
        // for all remaining SVs, including those which cross chromosomes,
        // assign them to the most precise cluster (ie that with the smallest positional range)
        for(SvClusterData currentVar : mUnassignedVariants) {

            for(int i = 0; i < 2; ++i) {

                boolean isStart = (i == 0);
                String isStartStr = isStart ? "start" : "end";

                // first start position, the end
                boolean found = false;
                for (SvCluster cluster : mClusters) {
                    SvCluster matchedCluster = cluster.findClosestCluster(currentVar.chromosome(isStart), currentVar.position(isStart));
                    if (matchedCluster != null) {
                        LOGGER.debug("variant: id({}) {}({}:{}) matched with cluster({}) on spanningVariant({})",
                                currentVar.id(), isStartStr, currentVar.chromosome(isStart), currentVar.position(isStart),
                                matchedCluster.getClusterId(), matchedCluster.getSpanningSV().posId());

                        if(isStart)
                            currentVar.setStartCluster(matchedCluster);
                        else
                            currentVar.setEndCluster(matchedCluster);

                        found = true;
                        break;
                    }
                }

                if (!found) {
                    LOGGER.debug("variant: id({}) {}({}:{}) unmatched",
                            currentVar.id(), isStartStr, currentVar.chromosome(isStart), currentVar.position(isStart));
                }
            }
        }
    }

    public void runClusteringByBaseDistance()
    {
        // purely for logging and verification
        logAllLinkedSVs();

        int currentIndex = 0;

        while(currentIndex < mUnassignedVariants.size())
        {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            // make a new cluster
            SvCluster newCluster = new SvCluster(getNextClusterId(), mClusteringUtils);

            LOGGER.debug("creating new cluster({}) with next variant({}:{})", newCluster.getClusterId(), currentIndex, currentVar.id());

            // first remove the current SV from consideration
            newCluster.addVariant(currentVar);
            mUnassignedVariants.remove(currentIndex); // index will remain the same and so point to the next item

            // and then search for all other linked ones
            findLinkedSVs(newCluster);

            mClusters.add(newCluster);

            LOGGER.debug("clusterCount({}) remainingVariants({})", mClusters.size(), mUnassignedVariants.size());
        }

        if(mConfig.getOutputCsvFile() != "")
            writeOutputCsvFile();
    }

    private int getNextClusterId() { return mNextClusterId++; }

    private void findLinkedSVs(SvCluster cluster) {
        // look for any other SVs which form part of this cluster
        int currentIndex = 0;

        while (currentIndex < mUnassignedVariants.size()) {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            // compare with all other SVs in this cluster
            boolean matched = false;
            for (SvClusterData otherVar : cluster.getSVs())
            {
                // test each possible linkage
                if (!mClusteringUtils.areVariantsLinked(currentVar, otherVar))
                {
                    LOGGER.debug("non-linked SVs: v1({}) and v2({})",
                            currentVar.posId(), otherVar.posId());
                    continue;
                }

                cluster.addVariant(currentVar);
                LOGGER.debug("cluster({}) add matched variant({}), totalCount({})", cluster.getClusterId(), currentVar.id(), cluster.getSVs().size());

                matched = true;
                break;
            }

            if(matched)
            {
                mUnassignedVariants.remove(currentIndex);

                // as soon as a new SV is added to this cluster, need to start checking from the beginning again
                currentIndex = 0;
            }
            else
            {
                ++currentIndex;
            }
        }
    }

    private void logAllLinkedSVs()
    {
        // purely for validation purposes
        for (final SvClusterData v1 : mAllVariants) {
            for (final SvClusterData v2 : mAllVariants) {

                if(mClusteringUtils.areVariantsLinked(v1, v2))
                {
                    LOGGER.debug("linked SVs: v1({}) and v2({})",
                            v1.posId(), v2.posId());
                }
            }
        }
    }

    private void writeOutputCsvFile()
    {
        String outputFileName = mConfig.getOutputCsvFile();

        Path outputFile = Paths.get(outputFileName + ".csv");

        try (final BufferedWriter writer = Files.newBufferedWriter(outputFile)) {

            writer.write("SAMPLE,CLUSTER_ID,FOOTPRINT_ID,SV_ID,TYPE,CRMS_START,POS_START,ORIENT_START,CRMS_END,POS_END,ORIENT_END\n");

            for(final SvCluster cluster : mClusters)
            {
                for(SvFootprint footprint : cluster.getFootprints()) {
                    for (final SvClusterData variant : footprint.getSVs()) {
                        writer.write(
                                String.format("%s,%d,%d,%s,%s,",
                                        mSampleId, cluster.getClusterId(), footprint.getFootprintId(), variant.id(), variant.type()));

                        writer.write(
                                String.format("%s,%d,%d,%s,%d,%d",
                                        variant.chromosome(true), variant.position(true), variant.orientation(true),
                                        variant.chromosome(false), variant.position(false), variant.orientation(false)
                                ));

                        writer.newLine();
                    }
                }
            }
            writer.close();
        }
        catch (final IOException e) {
            LOGGER.error("error writing outputFile{}", outputFileName);
        }
    }

    public List<SvCluster> getClusters() { return mClusters; }

    private final SvClusterData getSvData(final EnrichedStructuralVariant variant)
    {
        for(final SvClusterData svData : mAllVariants)
        {
            if(svData.id() == variant.id())
                return svData;
        }

        return null;
    }
}

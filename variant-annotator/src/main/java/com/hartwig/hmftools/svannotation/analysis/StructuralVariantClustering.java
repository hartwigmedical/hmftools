package com.hartwig.hmftools.svannotation.analysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svannotation.FragileSiteAnnotator;
import com.hartwig.hmftools.svannotation.LineElementAnnotator;
import com.hartwig.hmftools.svannotation.SvPONAnnotator;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.List;

public class StructuralVariantClustering {

    private final SvClusteringConfig mConfig;
    private final SvUtilities mClusteringUtils;
    private final ClusterAnalyser mAnalyser;

    // data per run (ie sample)
    private String mSampleId;
    private List<SvClusterData> mAllVariants; // the original list to analyse
    private List<SvClusterData> mUnassignedVariants; // the working list to be assigned
    private int mNextClusterId;

    private List<SvCluster> mClusters;

    BufferedWriter mFileWriter;
    SvPONAnnotator mSvPONAnnotator;
    FragileSiteAnnotator mFragileSiteAnnotator;
    LineElementAnnotator mLineElementAnnotator;

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantAnalyzer.class);

    public StructuralVariantClustering(final SvClusteringConfig config)
    {
        mConfig = config;
        mClusteringUtils = new SvUtilities(mConfig.getClusterBaseDistance());
        mFileWriter = null;
        mAnalyser = new ClusterAnalyser(config, mClusteringUtils);

        mSvPONAnnotator = new SvPONAnnotator();
        mSvPONAnnotator.loadPonFile(mConfig.getSvPONFile());

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.getFragileSiteFile());

        mLineElementAnnotator = new LineElementAnnotator();
        mLineElementAnnotator.loadLineElementsFile(mConfig.getLineElementFile());

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

//        if(LOGGER.isDebugEnabled())
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

        setChromosomalArms();

        addExternalAnnotations();

        // for now exclude Inserts
        removeInserts();

        clusterByBaseDistance();

        for(SvCluster cluster : mClusters)
        {
            mAnalyser.analyserCluster(cluster);
            cluster.setUniqueBreakends();
        }

        //runClusteringByLocals();

        if(mConfig.getOutputCsvPath() != "")
            writeBaseDistanceOutput();
    }

    private void setChromosomalArms()
    {
        for (SvClusterData var : mAllVariants)
        {
            String startArm = mClusteringUtils.getChromosomalArm(var.chromosome(true), var.position(true));
            String endArm = mClusteringUtils.getChromosomalArm(var.chromosome(false), var.position(false));
            var.setChromosomalArms(startArm, endArm);
        }
    }

    private void addExternalAnnotations()
    {
        for (SvClusterData var : mAllVariants)
        {
            mSvPONAnnotator.setPonOccurenceCount(var);
            var.setFragileSites(mFragileSiteAnnotator.isFragileSite(var, true), mFragileSiteAnnotator.isFragileSite(var, false));
            var.setLineElements(mLineElementAnnotator.isLineElement(var, true), mLineElementAnnotator.isLineElement(var, false));
        }
    }


    private void removeInserts()
    {
        int currentIndex = 0;

        while(currentIndex < mUnassignedVariants.size()) {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            // for now exclude Inserts
            if (currentVar.type() == StructuralVariantType.INS) {
                mUnassignedVariants.remove(currentIndex);
                continue;
            }

            ++currentIndex;
        }
    }

    public void clusterByBaseDistance()
    {
        // purely for logging and verification
        // logAllLinkedSVs();

        int currentIndex = 0;

        while(currentIndex < mUnassignedVariants.size()) {
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
                    //LOGGER.debug("non-linked SVs: v1({}) and v2({})", currentVar.posId(), otherVar.posId());
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

    private void writeBaseDistanceOutput()
    {
        try {

            BufferedWriter writer = null;

            if(mConfig.getUseCombinedOutputFile() && mFileWriter != null)
            {
                // check if can continue appending to an existing file
                writer = mFileWriter;
            }
            else
            {
                String outputFileName = mConfig.getOutputCsvPath();

                if(!outputFileName.endsWith("/"))
                    outputFileName += "/";

                if(mConfig.getUseCombinedOutputFile())
                    outputFileName += "CLUSTER.csv";
                else
                    outputFileName += mSampleId + ".csv";

                Path outputFile = Paths.get(outputFileName);

                boolean fileExists = Files.exists(outputFile);

                writer = Files.newBufferedWriter(outputFile, fileExists ? StandardOpenOption.APPEND : StandardOpenOption.CREATE_NEW);
                mFileWriter = writer;

                if(!fileExists) {
                    writer.write("SampleId,ClusterId,ClusterCount,Id,Type,Ploidy,PONCount,PONRegionCount,");
                    writer.write("ChrStart,PosStart,OrientStart,ArmStart,AdjAFStart,AdjCNStart,AdjCNChgStart,");
                    writer.write("ChrEnd,PosEnd,OrientEnd,ArmEnd,AdjAFEnd,AdjCNEnd,AdjCNChgEnd,");
                    writer.write("FSStart,FSEnd,FSCount,LEStart,LEEnd,LECount,DupBEStart,DupBEEnd,DupBECount,DupBESiteCount,");
                    writer.write("Desc,Consistency,ArmCount,IsTI,TICount,TILens,IsDSB,DSBCount,DSBLens,Annotations");
                    writer.write("\n");
                }
            }

            for(final SvCluster cluster : mClusters)
            {
                int lineElementCount = cluster.getLineElementCount();
                int fragileSiteCount = cluster.getFragileSiteCount();
                int duplicateBECount = cluster.getDuplicateBECount();
                int duplicateBESiteCount = cluster.getDuplicateBESiteCount();

                for (final SvClusterData var : cluster.getSVs()) {
                    writer.write(
                            String.format("%s,%d,%d,%s,%s,%.2f,%d,%d,",
                                    mSampleId, cluster.getClusterId(), cluster.getSVs().size(), var.id(), var.type(),
                                    var.getSvData().ploidy(), var.getPonCount(), var.getPonRegionCount()));

                    writer.write(
                            String.format("%s,%d,%d,%s,%.2f,%.2f,%.2f,%s,%d,%d,%s,%.2f,%.2f,%.2f,",
                                    var.chromosome(true), var.position(true), var.orientation(true), var.getStartArm(),
                                    var.getSvData().adjustedStartAF(), var.getSvData().adjustedStartCopyNumber(), var.getSvData().adjustedStartCopyNumberChange(),
                                    var.chromosome(false), var.position(false), var.orientation(false), var.getEndArm(),
                                    var.getSvData().adjustedEndAF(), var.getSvData().adjustedEndCopyNumber(), var.getSvData().adjustedEndCopyNumberChange()
                                    ));

                    writer.write(
                            String.format("%s,%s,%d,%s,%s,%d,%s,%s,%d,%d,",
                                    var.isStartFragileSite(), var.isEndFragileSite(), fragileSiteCount, var.isStartLineElement(), var.isEndLineElement(), lineElementCount,
                                    cluster.isSvDuplicateBE(var, true), cluster.isSvDuplicateBE(var, false), duplicateBECount, duplicateBESiteCount
                                    ));

                    writer.write(
                            String.format("%s,%d,%d,%s,%d,%s,%s,%d,%s,%s",
                                    cluster.getDesc(), cluster.getConsistencyCount(), cluster.getChromosomalArmCount(),
                                    cluster.isReplicationEvent(), cluster.getTICount(), cluster.getTempInsertLengths(),
                                    cluster.isDSBEvent(), cluster.getDSBCount(), cluster.getDSBLengths(), cluster.getAnnotations()
                                    ));

                    writer.newLine();
                }
            }
            // writer.close();

        }
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
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

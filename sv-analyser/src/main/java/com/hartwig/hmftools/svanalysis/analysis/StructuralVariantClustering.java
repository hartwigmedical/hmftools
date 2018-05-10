package com.hartwig.hmftools.svanalysis.analysis;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.annotators.ExternalSVAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

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
    ExternalSVAnnotator mExternalAnnotator;
    SvClusteringMethods mClusteringMethods;

    PerformanceCounter mPerfCounter;
    PerformanceCounter mPC1;
    PerformanceCounter mPC2;
    PerformanceCounter mPC3;
    PerformanceCounter mPC4;
    PerformanceCounter mPC5;
    PerformanceCounter mPC6;

    private static final Logger LOGGER = LogManager.getLogger(StructuralVariantClustering.class);

    public StructuralVariantClustering(final SvClusteringConfig config)
    {
        mConfig = config;
        mClusteringUtils = new SvUtilities(mConfig.getClusterBaseDistance());
        mFileWriter = null;
        mAnalyser = new ClusterAnalyser(config, mClusteringUtils);

        mClusteringMethods = new SvClusteringMethods(mClusteringUtils);

        mSvPONAnnotator = new SvPONAnnotator();
        mSvPONAnnotator.loadPonFile(mConfig.getSvPONFile());

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.getFragileSiteFile());

        mLineElementAnnotator = new LineElementAnnotator();
        mLineElementAnnotator.loadLineElementsFile(mConfig.getLineElementFile());

        mExternalAnnotator = new ExternalSVAnnotator();
        mExternalAnnotator.loadFile(mConfig.getExternalAnnotationsFile());

        mPerfCounter = new PerformanceCounter("Total");
        mPC1 = new PerformanceCounter("Annotate&Filter");
        // mPC2 = new PerformanceCounter("ExtAnn");
        mPC3 = new PerformanceCounter("BaseDist");
        mPC4 = new PerformanceCounter("Analyse");
        mPC5 = new PerformanceCounter("Nearest");
        mPC6 = new PerformanceCounter("CSV");

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
            mAllVariants.add(SvClusterData.from(enrichedSV));
        }
    }

    public void loadFromDatabase(final String sampleId, final List<SvClusterData> variants)
    {
        if (variants.isEmpty())
            return;

        clearState();

        mSampleId = sampleId;
        mAllVariants = Lists.newArrayList(variants);

        LOGGER.debug("loaded {} SVs", mAllVariants.size());

//        if(LOGGER.isDebugEnabled())
//        for(SvClusterData sv : variants) {
//            LOGGER.debug("adding SV: id({}) start({}:{}) end({}:{})",
//                    sv.id(), sv.chromosome(true), sv.position(true), sv.chromosome(false), sv.position(false));
//        }
    }

    public void runClustering()
    {
        if(mAllVariants.isEmpty())
            return;

        mPerfCounter.start();

        mPC1.start();

        annotateAndFilterVariants();

        mPC1.stop();

//        mPC2.start();
//        // LOGGER.debug("sample({}) clustering {} variants", mSampleId, mAllVariants.size());
//        mPC2.stop();

        mPC3.start();

        LOGGER.debug("sample({}) clustering {} variants", mSampleId, mAllVariants.size());

        clusterByBaseDistance();

        mPC3.stop();
        mPC4.start();

        for(SvCluster cluster : mClusters)
        {
            mAnalyser.setClusterStats(cluster);
            cluster.setUniqueBreakends();
            mAnalyser.findLinkedPairs(mSampleId, cluster);
            mAnalyser.findSvChains(mSampleId, cluster);
            mAnalyser.resolveTransitiveSVs(mSampleId, cluster);
        }

        mPC4.stop();
        mPC5.start();

//        if(!mExternalAnnotator.hasExternalData()) {
//            // skip this step if sourced externally
//            mClusteringMethods.setClosestSVData(mAllVariants, mClusters);
//            mClusteringMethods.setClosestLinkedSVData(mAllVariants);
//        }
        mClusteringMethods.setChromosomalArmStats(mAllVariants);

        mPC5.stop();
        mPC6.start();

        if(mConfig.getOutputCsvPath() != "")
            writeClusterDataOutput();

        mPC6.stop();
        mPerfCounter.stop();
    }

    private void annotateAndFilterVariants()
    {
        int currentIndex = 0;

        while(currentIndex < mAllVariants.size()) {

            SvClusterData var = mAllVariants.get(currentIndex);

            // for now exclude Inserts
            if (var.type() == StructuralVariantType.INS) {
                mAllVariants.remove(currentIndex);
                continue;
            }

            mSvPONAnnotator.setPonOccurenceCount(var);

            // exclude PON
            if (var.getPonCount() >= 2) {
                mAllVariants.remove(currentIndex);
                continue;
            }

            var.setFragileSites(mFragileSiteAnnotator.isFragileSite(var, true), mFragileSiteAnnotator.isFragileSite(var, false));
            var.setLineElements(mLineElementAnnotator.isLineElement(var, true), mLineElementAnnotator.isLineElement(var, false));

            String startArm = mClusteringUtils.getChromosomalArm(var.chromosome(true), var.position(true));
            String endArm = mClusteringUtils.getChromosomalArm(var.chromosome(false), var.position(false));
            var.setChromosomalArms(startArm, endArm);

            ++currentIndex;
        }
    }

    public void clusterByBaseDistance()
    {
        // purely for logging and verification
        // logAllLinkedSVs();

        mUnassignedVariants = Lists.newArrayList(mAllVariants);

        // assign each variant once to a cluster using proximity as a test
        int currentIndex = 0;

        while(currentIndex < mUnassignedVariants.size()) {
            SvClusterData currentVar = mUnassignedVariants.get(currentIndex);

            // make a new cluster
            SvCluster newCluster = new SvCluster(getNextClusterId(), mClusteringUtils);

            LOGGER.debug("creating new cluster({}) with next variant({}:{})", newCluster.getId(), currentIndex, currentVar.id());

            // first remove the current SV from consideration
            newCluster.addVariant(currentVar);
            mUnassignedVariants.remove(currentIndex); // index will remain the same and so point to the next item

            // and then search for all other linked ones
            mClusteringMethods.findLinkedSVsByDistance(newCluster, mUnassignedVariants);

            mClusters.add(newCluster);

            LOGGER.debug("clusterCount({}) remainingVariants({})", mClusters.size(), mUnassignedVariants.size());
        }
    }

    private int getNextClusterId() { return mNextClusterId++; }

    private void logAllLinkedSVs()
    {
        // purely for validation purposes
        for (final SvClusterData v1 : mAllVariants) {
            for (final SvClusterData v2 : mAllVariants) {

                if(mClusteringUtils.areVariantsLinkedByDistance(v1, v2))
                {
                    LOGGER.debug("linked SVs: v1({}) and v2({})",
                            v1.posId(), v2.posId());
                }
            }
        }
    }

    private void writeClusterDataOutput()
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
                    // SV info
                    writer.write("SampleId,ClusterId,ClusterCount,Id,Type,Ploidy,PONCount,PONRegionCount,");
                    writer.write("ChrStart,PosStart,OrientStart,ArmStart,AdjAFStart,AdjCNStart,AdjCNChgStart,");
                    writer.write("ChrEnd,PosEnd,OrientEnd,ArmEnd,AdjAFEnd,AdjCNEnd,AdjCNChgEnd,Homology,");
                    writer.write("FSStart,FSEnd,LEStart,LEEnd,DupBEStart,DupBEEnd,");
                    writer.write("ArmCountStart,ArmExpStart,ArmCountEnd,ArmExpEnd,");

                    // cluster-level info
                    writer.write("FSCount,LECount,DupBECount,DupBESiteCount,");
                    writer.write("ClusterDesc,Consistency,ArmCount,");

                    // linked pair info
                    writer.write("LnkStartSv,LnkStartType,LnkStartLen,LnkEndSv,LnkEndType,LnkEndLen,");

                    // chain info
                    writer.write("ChainId,ChainCount,ChainTICount,ChainDBCount,ChainIndex,");

                    // transitive info
                    writer.write("TransType,TransLen,TransSvLinks,");

                    // external annotations if any
                    if(mExternalAnnotator.hasExternalData())
                        writer.write(String.format("%s", mExternalAnnotator.getFieldNames()));
                    else
                        writer.write("NearestLen,NearestType,NearestTILen,NearestDBLen");

                    writer.newLine();
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
                                    mSampleId, cluster.getId(), cluster.getCount(), var.id(), var.type(),
                                    var.getSvData().ploidy(), var.getPonCount(), var.getPonRegionCount()));

                    writer.write(
                            String.format("%s,%d,%d,%s,%.2f,%.2f,%.2f,%s,%d,%d,%s,%.2f,%.2f,%.2f,%s,",
                                    var.chromosome(true), var.position(true), var.orientation(true), var.getStartArm(),
                                    var.getSvData().adjustedStartAF(), var.getSvData().adjustedStartCopyNumber(), var.getSvData().adjustedStartCopyNumberChange(),
                                    var.chromosome(false), var.position(false), var.orientation(false), var.getEndArm(),
                                    var.getSvData().adjustedEndAF(), var.getSvData().adjustedEndCopyNumber(), var.getSvData().adjustedEndCopyNumberChange(),
                                    var.getSvData().insertSequence().isEmpty() ? var.getSvData().homology() : ""
                                    ));

                    writer.write(
                            String.format("%s,%s,%s,%s,%s,%s,%s,",
                                    // var.getNearestSVLength(), var.getNearestSVLinkType(), var.getNearestTILength(), var.getNearestDBLength(),
                                    var.isStartFragileSite(), var.isEndFragileSite(),
                                    var.isStartLineElement(), var.isEndLineElement(),
                                    var.isDupBEStart(), var.isDupBEEnd(),
                                    mClusteringMethods.getChrArmData(var)));

                    writer.write(
                            String.format("%d,%d,%d,%d,",
                                    fragileSiteCount, lineElementCount, duplicateBECount, duplicateBESiteCount));

                    writer.write(
                            String.format("%s,%d,%d,",
                                    cluster.getDesc(), cluster.getConsistencyCount(), cluster.getChromosomalArmCount()));

                    // linked pair info
                    SvLinkedPair startLP = cluster.findLinkedPair(var, true);
                    if(startLP != null)
                        writer.write(String.format("%s,%s,%d,",
                                startLP.first().equals(var) ? startLP.second().id() : startLP.first().id(),
                                startLP.linkType(), startLP.length()));
                    else
                        writer.write("0,,-1,");

                    SvLinkedPair endLP = cluster.findLinkedPair(var, false);
                    if(endLP != null)
                        writer.write(String.format("%s,%s,%d,",
                                endLP.first().equals(var) ? endLP.second().id() : endLP.first().id(),
                                endLP.linkType(), endLP.length()));
                    else
                        writer.write("0,,-1,");

                    // chain info
                    SvChain chain = cluster.findChain(var);

                    if(chain != null)
                        writer.write(String.format("%d,%d,%d,%d,%d,",
                                chain.getId(), chain.getLinkCount(), chain.getTICount(), chain.getDBCount(), chain.getSvIndex(var)));
                    else
                        writer.write("0,0,0,0,0,");

                    writer.write(String.format("%s,%d,%s,", var.getTransType(), var.getTransLength(), var.getTransSvLinks()));

                    if(mExternalAnnotator.hasExternalData()) {
                        writer.write(String.format("%s", mExternalAnnotator.getSVData(var)));
                    }
                    else
                    {
                        writer.write(String.format("%d,%s,%d,%d",
                                var.getNearestSVLength(), var.getNearestSVLinkType(),
                                var.getNearestTILength(), var.getNearestDBLength()));
                    }

                    writer.newLine();
                }
            }

        }
        catch (final IOException e) {
            LOGGER.error("error writing to outputFile");
        }
    }

    public void close()
    {
        if(mFileWriter == null)
            return;

        try
        {
            mFileWriter.close();
        }
        catch (final IOException e) {
        }

        // log perf stats
        mPerfCounter.logStats();
        mPC1.logStats();
        //mPC2.logStats();
        mPC3.logStats();
        mPC4.logStats();
        mPC5.logStats();
        mPC6.logStats();
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

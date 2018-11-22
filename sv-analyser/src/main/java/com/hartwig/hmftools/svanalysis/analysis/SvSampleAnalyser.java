package com.hartwig.hmftools.svanalysis.analysis;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator.REGION_DISTANCE;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_LOW_QUALITY;
import static com.hartwig.hmftools.svanalysis.types.SvCluster.RESOLVED_TYPE_SIMPLE_SV;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_ASMB_ONLY;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.annotators.ExternalSVAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator;
import com.hartwig.hmftools.common.variant.structural.annotation.SvPONAnnotator;
import com.hartwig.hmftools.svanalysis.types.SvArmGroup;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvCluster;
import com.hartwig.hmftools.common.variant.structural.annotation.SvPON;
import com.hartwig.hmftools.svanalysis.types.SvVarData;
import com.hartwig.hmftools.svanalysis.types.SvLinkedPair;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SvSampleAnalyser {

    private final SvClusteringConfig mConfig;
    private final SvUtilities mClusteringUtils;
    private final ClusterAnalyser mAnalyser;

    // data per run (ie sample)
    private String mSampleId;
    private List<SvVarData> mAllVariants; // the original list to analyse

    BufferedWriter mFileWriter;
    SvPONAnnotator mSvPONAnnotator;
    FragileSiteAnnotator mFragileSiteAnnotator;
    LineElementAnnotator mLineElementAnnotator;
    ExternalSVAnnotator mExternalAnnotator;
    SvClusteringMethods mClusteringMethods;
    CNAnalyser mCnAnalyser;

    PerformanceCounter mPerfCounter;
    PerformanceCounter mPc1;
    PerformanceCounter mPc2;
    PerformanceCounter mPc3;
    PerformanceCounter mPc4;
    PerformanceCounter mPc5;

    private static final Logger LOGGER = LogManager.getLogger(SvSampleAnalyser.class);

    public SvSampleAnalyser(final SvClusteringConfig config)
    {
        mConfig = config;
        mClusteringUtils = new SvUtilities(mConfig.ProximityDistance);
        mClusteringMethods = new SvClusteringMethods(mClusteringUtils);
        mAnalyser = new ClusterAnalyser(config, mClusteringUtils, mClusteringMethods);

        mFileWriter = null;

        mSvPONAnnotator = new SvPONAnnotator();
        mSvPONAnnotator.loadPonFile(mConfig.SvPONFile);

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.FragileSiteFile);

        mLineElementAnnotator = new LineElementAnnotator();
        mLineElementAnnotator.loadLineElementsFile(mConfig.LineElementFile);

        mExternalAnnotator = new ExternalSVAnnotator();
        mExternalAnnotator.loadFile(mConfig.ExternalAnnotationsFile);

        mCnAnalyser = null;
        if(!mConfig.LOHDataFile.isEmpty())
        {
            mCnAnalyser = new CNAnalyser(mConfig.OutputCsvPath, null);
            mCnAnalyser.loadLOHFromCSV(mConfig.LOHDataFile, "");
            mClusteringMethods.setSampleLohData(mCnAnalyser.getSampleLohData());
        }

        mPerfCounter = new PerformanceCounter("Total");

        mPc1 = new PerformanceCounter("Annotate&Filter");
        mPc2 = new PerformanceCounter("ArmsStats");
        mPc3 = new PerformanceCounter("ClusterAndAnalyse");
        // mPc4 = new PerformanceCounter("Analyse");
        mPc5 = new PerformanceCounter("WriteCSV");

        mPerfCounter.start();

        clearState();
    }

    public final List<SvCluster> getClusters() { return mAnalyser.getClusters(); }

    private void clearState()
    {
        mSampleId = "";
        mAllVariants = Lists.newArrayList();
    }

    public void loadFromDatabase(final String sampleId, final List<SvVarData> variants)
    {
        clearState();

        if (variants.isEmpty())
            return;

        mSampleId = sampleId;
        mAllVariants = Lists.newArrayList(variants);

        LOGGER.debug("loaded {} SVs", mAllVariants.size());
    }

    public void analyse()
    {
        if(mAllVariants.isEmpty())
            return;

        mPc1.start();
        annotateAndFilterVariants();
        mPc1.stop();

        LOGGER.debug("sample({}) clustering {} variants", mSampleId, mAllVariants.size());

        mPc2.start();
        mClusteringMethods.setChromosomalArmStats(mAllVariants);
        mClusteringMethods.populateChromosomeBreakendMap(mAllVariants);
        // mClusteringMethods.calcCopyNumberData(mSampleId);
        // mClusteringMethods.createCopyNumberSegments();
        mClusteringMethods.annotateNearestSvData();
        LinkFinder.findDeletionBridges(mClusteringMethods.getChrBreakendMap());
        mClusteringMethods.setSimpleVariantLengths(mSampleId);
        mLineElementAnnotator.setSuspectedLineElements(mClusteringMethods.getChrBreakendMap(), mConfig.ProximityDistance);
        mPc2.stop();

        mPc3.start();

        mAnalyser.setSampleData(mSampleId, mAllVariants);
        mAnalyser.clusterAndAnalyse();

        mPc3.stop();

        /*
        mClusteringMethods.clusterByBaseDistance(mAllVariants, mClusters);
        mAnalyser.setSampleId(mSampleId);
        mAnalyser.findSimpleCompleteChains();
        mClusteringMethods.mergeClusters(mSampleId, mClusters);

        // log basic clustering details
        for(SvCluster cluster : mClusters)
        {
            if(cluster.getCount() > 1)
                cluster.logDetails();
        }

        mPc4.start();
        mAnalyser.findLinksAndChains();
        // mClusteringMethods.logInversionPairData(mSampleId, mClusters);
        mPc4.stop();
        */

        logSampleClusterInfo();

        mPc5.start();

        if(!mConfig.OutputCsvPath.isEmpty())
            writeClusterDataOutput();

        mPc5.stop();
    }

    private void annotateAndFilterVariants()
    {
        int currentIndex = 0;

        while(currentIndex < mAllVariants.size())
        {
            SvVarData var = mAllVariants.get(currentIndex);

            if(mExternalAnnotator.hasData())
            {
                mExternalAnnotator.setSVData(mSampleId, var);
            }
            else
            {
                setPonOccurenceCount(var);
                var.setFragileSites(mFragileSiteAnnotator.isFragileSite(var, true), mFragileSiteAnnotator.isFragileSite(var, false));
                var.setLineElement(mLineElementAnnotator.isLineElement(var, true), true);
                var.setLineElement(mLineElementAnnotator.isLineElement(var, false), false);
            }

            // exclude PON
            if (var.getPonCount() >= 2)
            {
                LOGGER.info("filtering sv({}) with PonCount({})", var.id(), var.getPonCount());
                mAllVariants.remove(currentIndex);
                continue;
            }

            String startArm = mClusteringUtils.getChromosomalArm(var.chromosome(true), var.position(true));

            String endArm = "";
            if(!var.isNullBreakend())
                endArm = mClusteringUtils.getChromosomalArm(var.chromosome(false), var.position(false));
            else
                endArm = CHROMOSOME_ARM_P;

            var.setChromosomalArms(startArm, endArm);

            ++currentIndex;
        }
    }

    private void writeClusterDataOutput()
    {
        try
        {
            BufferedWriter writer = null;

            if(mConfig.UseCombinedOutputFile && mFileWriter != null)
            {
                // check if can continue appending to an existing file
                writer = mFileWriter;
            }
            else
            {
                String outputFileName = mConfig.OutputCsvPath;

                if(!outputFileName.endsWith("/"))
                    outputFileName += "/";

                if(mConfig.UseCombinedOutputFile)
                    outputFileName += "CLUSTER.csv";
                else
                    outputFileName += mSampleId + ".csv";

                Path outputFile = Paths.get(outputFileName);

                writer = Files.newBufferedWriter(outputFile, StandardOpenOption.CREATE);
                mFileWriter = writer;

                // definitional fields
                writer.write("SampleId,Id,Type,ChrStart,PosStart,OrientStart,ChrEnd,PosEnd,OrientEnd");

                // position and copy number
                writer.write(",ArmStart,AdjAFStart,AdjCNStart,AdjCNChgStart,ArmEnd,AdjAFEnd,AdjCNEnd,AdjCNChgEnd,Ploidy");

                // cluster info
                writer.write(",ClusterId,SubClusterId,ClusterCount,ClusterReason");

                // cluster-level info
                writer.write(",ClusterDesc,IsResolved,ResolvedType,Consistency,ArmCount");

                // SV info
                writer.write(",Homology,InexactHOStart,InexactHOEnd,InsertSeq,Imprecise,PONCount,QualScore");

                // location attributes
                writer.write(",FSStart,FSEnd,LEStart,LEEnd,DupBEStart,DupBEEnd,ArmCountStart,ArmExpStart,ArmCountEnd,ArmExpEnd");

                // linked pair info
                writer.write(",LnkSvStart,LnkTypeStart,LnkLenStart,LnkSvEnd,LnkTypeEnd,LnkLenEnd");

                // GRIDDS caller info
                writer.write(",AsmbStart,AsmbEnd,AsmbMatchStart,AsmbMatchEnd");

                // chain info
                writer.write(",ChainId,ChainCount,ChainIndex");

                // proximity info and other link info
                writer.write(",NearestLen,NearestType,DBLenStart,DBLenEnd,FbkLnkStart,FbkLenStart,FbkLnkEnd,FbkLenEnd");

                // transitive info
                // writer.write(",TransType,TransLen,TransSvLinks");

                writer.newLine();
            }

            int lineCount = 0;
            int svCount = 0;

            for(final SvVarData var : mAllVariants)
            {
                final SvCluster cluster = var.getCluster();

                if(cluster == null)
                {
                    LOGGER.error("SV({}) not assigned to any cluster", var.posId());
                    continue;
                }

                int clusterSvCount = cluster.getUniqueSvCount();

                SvCluster subCluster = cluster;
                if(cluster.hasSubClusters())
                {
                    for(final SvCluster sc : cluster.getSubClusters())
                    {
                        if(sc.getSVs().contains(var))
                        {
                            subCluster = sc;
                            break;
                        }
                    }
                }

                final StructuralVariantData dbData = var.getSvData();

                ++svCount;

                writer.write(
                        String.format("%s,%s,%s,%s,%d,%d,%s,%d,%d",
                                mSampleId, var.id(), var.typeStr(),
                                var.chromosome(true), var.position(true), var.orientation(true),
                                var.chromosome(false), var.position(false), var.orientation(false)));

                writer.write(
                        String.format(",%s,%.2f,%.2f,%.2f,%s,%.2f,%.2f,%.2f,%.2f",
                                var.arm(true), dbData.adjustedStartAF(), dbData.adjustedStartCopyNumber(), dbData.adjustedStartCopyNumberChange(),
                                var.arm(false), dbData.adjustedEndAF(), dbData.adjustedEndCopyNumber(), dbData.adjustedEndCopyNumberChange(), dbData.ploidy()));

                writer.write(
                        String.format(",%d,%d,%d,%s",
                                cluster.getId(), subCluster.getId(), clusterSvCount, var.getClusterReason()));

                writer.write(
                        String.format(",%s,%s,%s,%d,%d",
                                cluster.getDesc(), cluster.isResolved(), cluster.getResolvedType(), cluster.getConsistencyCount(), cluster.getChromosomalArmCount()));

                writer.write(
                        String.format(",%s,%d,%d,%s,%s,%d,%.0f",
                                dbData.insertSequence().isEmpty() && var.type() != StructuralVariantType.INS ? dbData.homology() : "",
                                dbData.inexactHomologyOffsetStart(), dbData.inexactHomologyOffsetEnd(),
                                dbData.insertSequence(), dbData.imprecise(), var.getPonCount(), dbData.qualityScore()));

                writer.write(
                        String.format(",%s,%s,%s,%s,%s,%s,%s",
                                var.isFragileSite(true), var.isFragileSite(false),
                                var.getLineElement(true), var.getLineElement(false),
                                var.isDupBreakend(true), var.isDupBreakend(false),
                                mClusteringMethods.getChrArmData(var)));

                // linked pair info
                final SvLinkedPair startLP = var.getLinkedPair(true) != null ? var.getLinkedPair(true) : cluster.getLinkedPair(var, true);
                String startLinkStr = "0,,-1";
                String assemblyMatchStart = var.getAssemblyMatchType(true);
                if(startLP != null)
                {
                    startLinkStr = String.format("%s,%s,%d",
                            startLP.first().equals(var, true) ? startLP.second().origId() : startLP.first().origId(),
                            startLP.linkType(), startLP.length());
                }

                final SvLinkedPair endLP = var.getLinkedPair(false) != null ? var.getLinkedPair(false) : cluster.getLinkedPair(var, false);
                String endLinkStr = "0,,-1";
                String assemblyMatchEnd = var.getAssemblyMatchType(false);
                if(endLP != null)
                {
                    endLinkStr = String.format("%s,%s,%d",
                            endLP.first().equals(var, true) ? endLP.second().origId() : endLP.first().origId(),
                            endLP.linkType(), endLP.length());
                }

                if(assemblyMatchStart.equals(ASSEMBLY_MATCH_ASMB_ONLY) || assemblyMatchEnd.equals(ASSEMBLY_MATCH_ASMB_ONLY))
                {
                    LOGGER.debug("sample({}) var({}) has unlinked assembly TIs", mSampleId, var.posId());
                }

                // assembly info
                writer.write(String.format(",%s,%s,%s,%s,%s,%s",
                        startLinkStr, endLinkStr, var.getAssemblyData(true), var.getAssemblyData(false), assemblyMatchStart, assemblyMatchEnd));

                // chain info
                final SvChain chain = cluster.findChain(var);
                String chainStr = ",0,0,";

                if(chain != null)
                {
                    chainStr = String.format(",%d,%d,%s", chain.getId(), chain.getUniqueSvCount(), chain.getSvIndices(var));
                }

                writer.write(chainStr);

                writer.write(String.format(",%d,%s,%d,%d,%s,%d,%s,%d",
                        var.getNearestSvDistance(), var.getNearestSvRelation(),
                        var.getDBLink(true) != null ? var.getDBLink(true).length() : -1,
                        var.getDBLink(false) != null ? var.getDBLink(false).length() : -1,
                        var.getFoldbackLink(true), var.getFoldbackLen(true),
                        var.getFoldbackLink(false), var.getFoldbackLen(false)));

                // writer.write(String.format(",%s,%d,%s", var.getTransType(), var.getTransLength(), var.getTransSvLinks()));

                ++lineCount;
                writer.newLine();

                if(svCount != lineCount)
                {
                    LOGGER.error("inconsistent output");
                }
            }

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
        }
    }

    public void setPonOccurenceCount(SvVarData svData)
    {
        if(mSvPONAnnotator == null || !mSvPONAnnotator.hasEntries())
            return;

        if(mSvPONAnnotator.getPonList().isEmpty())
            return;

        final List<SvPON> ponList = mSvPONAnnotator.getPonList();
        final Map<String, Integer> chrIndexMap = mSvPONAnnotator.getChrIndexMap();

        // use CRMS to find start index
        if(!chrIndexMap.containsKey(svData.chromosome(true)))
            return;

        int index = chrIndexMap.get(svData.chromosome(true));

        for(; index < ponList.size(); ++index)
        {
            final SvPON svPon = ponList.get(index);

            if(!svPon.chrStart().equals(svData.chromosome(true)))
            {
                return;
            }

            if(isMatch(svData, svPon,0))
            {
                LOGGER.debug("var({}) found in PON with count({})", svData.posId(), svPon.count());
                svData.setPonCount(svPon.count());
                svData.setPonRegionCount(svPon.count()); // set to match
                return;
            }
            else if(svData.getPonRegionCount() == 0 && isMatch(svData, svPon, REGION_DISTANCE))
            {
                LOGGER.debug("var({}) found in PON region({}->{}) with count({})",
                        svData.posId(), svPon.posStart(), svPon.posEnd(), svPon.count());

                svData.setPonRegionCount(svPon.count());
            }
        }
    }

    private static boolean isMatch(final SvVarData svData, final SvPON svPon, final int permittedDiff)
    {
        if(!svPon.chrStart().equals(svData.chromosome(true))
        || !svPon.chrEnd().equals(svData.chromosome(false))
        || svPon.orientStart() != svData.orientation(true)
        || svPon.orientEnd() != svData.orientation(false)
        || !svPon.type().equals(svData.type().toString()))
        {
            return false;
        }

        return Math.abs(svPon.posStart() - svData.position(true)) <= permittedDiff
                && Math.abs(svPon.posEnd() - svData.position(false)) <= permittedDiff;
    }

    private void logSampleClusterInfo()
    {
        int simpleClusterCount = 0;
        int complexClusterCount = 0;
        Map<String, Integer> armSimpleClusterCount = new HashMap();
        Map<String, Integer> armComplexClusterCount = new HashMap();

        for(final SvCluster cluster : mAnalyser.getClusters())
        {
            Map<String, Integer> targetMap;

            if(cluster.getResolvedType() == RESOLVED_LOW_QUALITY)
                continue;

            if(cluster.getResolvedType() == RESOLVED_TYPE_SIMPLE_SV || cluster.isSyntheticSimpleType())
            {
                ++simpleClusterCount;
                targetMap = armSimpleClusterCount;
            }
            else
            {
                ++complexClusterCount;
                targetMap = armComplexClusterCount;
            }

            for(final SvArmGroup armGroup : cluster.getArmGroups())
            {
                if(targetMap.containsKey(armGroup))
                {
                    targetMap.put(armGroup.id(), targetMap.get(armGroup.id()) + 1);
                }
                else
                {
                    targetMap.put(armGroup.id(), 1);
                }
            }
        }

        int armsWithExcessComplexClusters = 0;

        for(final Map.Entry<String, Integer> entry : armComplexClusterCount.entrySet())
        {
            final String chrArm = entry.getKey();
            int complexCount = entry.getValue();
            Integer simpleCount = armSimpleClusterCount.containsKey(chrArm) ? armSimpleClusterCount.get(chrArm) : 0;

            if(simpleCount > 0 && complexCount > simpleCount)
            {
                LOGGER.debug("chrArm({}) clusters simple({}) vs complex({})", chrArm, simpleCount, complexCount);
                ++armsWithExcessComplexClusters;
            }
        }

        if(complexClusterCount > simpleClusterCount || armsWithExcessComplexClusters >= 2)
        {
            LOGGER.info("sample({}) clusters total({}) simple({}) complex({}) excessComplexArms({})",
                    mSampleId, mAnalyser.getClusters().size(), simpleClusterCount, complexClusterCount, armsWithExcessComplexClusters);
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
        catch (final IOException e)
        {
        }

        // log perf stats
        mPerfCounter.stop();
        mPerfCounter.logStats(false);
        mPc1.logStats(false);
        mPc2.logStats(false);
        mPc3.logStats(false);
        // mPc4.logStats(false);
        mPc5.logStats(false);

        mAnalyser.logStats();
    }

}

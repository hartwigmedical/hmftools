package com.hartwig.hmftools.svanalysis.analysis;

import static com.hartwig.hmftools.svanalysis.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_ASMB_ONLY;
import static com.hartwig.hmftools.svanalysis.types.SvLinkedPair.ASSEMBLY_MATCH_NONE;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;
import com.hartwig.hmftools.svanalysis.annotators.ExternalSVAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.FragileSiteAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.GeneAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.LineElementAnnotator;
import com.hartwig.hmftools.svanalysis.annotators.SvPONAnnotator;
import com.hartwig.hmftools.svanalysis.types.SvBreakend;
import com.hartwig.hmftools.svanalysis.types.SvChain;
import com.hartwig.hmftools.svanalysis.types.SvClusterData;
import com.hartwig.hmftools.svanalysis.types.SvGeneData;
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
    private List<SvClusterData> mAllVariants; // the original list to analyse

    private List<SvCluster> mClusters;

    BufferedWriter mFileWriter;
    SvPONAnnotator mSvPONAnnotator;
    FragileSiteAnnotator mFragileSiteAnnotator;
    LineElementAnnotator mLineElementAnnotator;
    ExternalSVAnnotator mExternalAnnotator;
    SvClusteringMethods mClusteringMethods;
    GeneAnnotator mGeneAnnotator;

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
        mClusteringUtils = new SvUtilities(mConfig.ClusterBaseDistance);
        mFileWriter = null;
        mAnalyser = new ClusterAnalyser(config, mClusteringUtils);

        mClusteringMethods = new SvClusteringMethods(mClusteringUtils);

        mSvPONAnnotator = new SvPONAnnotator();
        mSvPONAnnotator.loadPonFile(mConfig.SvPONFile);

        mFragileSiteAnnotator = new FragileSiteAnnotator();
        mFragileSiteAnnotator.loadFragileSitesFile(mConfig.FragileSiteFile);

        mLineElementAnnotator = new LineElementAnnotator();
        mLineElementAnnotator.loadLineElementsFile(mConfig.LineElementFile);

        mExternalAnnotator = new ExternalSVAnnotator();
        mExternalAnnotator.loadFile(mConfig.ExternalAnnotationsFile);

        mGeneAnnotator = new GeneAnnotator();
        mGeneAnnotator.loadGeneDriverFile(mConfig.GeneDataFile);

        mPerfCounter = new PerformanceCounter("Total");

        mPc1 = new PerformanceCounter("Annotate&Filter");
        mPc2 = new PerformanceCounter("ArmsStats");
        mPc3 = new PerformanceCounter("Clusters");
        mPc4 = new PerformanceCounter("Analyse");
        mPc5 = new PerformanceCounter("WriteCSV");

        mPerfCounter.start();

        clearState();
    }

    private void clearState()
    {
        mSampleId = "";
        mAllVariants = Lists.newArrayList();
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
        mClusteringMethods.calcCopyNumberData(mSampleId);
        mClusteringMethods.annotateNearestSvData();
        mPc2.stop();

        mPc3.start();
        mClusteringMethods.clusterByBaseDistance(mAllVariants, mClusters);
        mClusteringMethods.mergeClusters(mClusters);
        mPc3.stop();

        mPc4.start();

        for(SvCluster cluster : mClusters)
        {
            mAnalyser.setClusterStats(cluster);
            cluster.setUniqueBreakends();
        }

        mAnalyser.setChrCopyNumberData(mClusteringMethods.getChrCopyNumberMap());
        mAnalyser.setChrBreakendData(mClusteringMethods.getChrBreakendMap());
        mAnalyser.findLinksAndChains(mSampleId, mClusters);

        mPc4.stop();

        if(mGeneAnnotator.hasData())
        {
            for (SvClusterData var : mAllVariants)
            {
                mGeneAnnotator.addGeneData(mSampleId, var);
            }

            mGeneAnnotator.reportGeneMatchData(mSampleId);
        }

        mPc5.start();

        if(!mConfig.OutputCsvPath.isEmpty())
            writeClusterDataOutput();

        mPc5.stop();
    }

    private void annotateAndFilterVariants()
    {
        int currentIndex = 0;

        while(currentIndex < mAllVariants.size()) {

            SvClusterData var = mAllVariants.get(currentIndex);

            if(mExternalAnnotator.hasData())
            {
                mExternalAnnotator.setSVData(mSampleId, var);
            }
            else {

                mSvPONAnnotator.setPonOccurenceCount(var);
                var.setFragileSites(mFragileSiteAnnotator.isFragileSite(var, true), mFragileSiteAnnotator.isFragileSite(var, false));
                var.setLineElements(mLineElementAnnotator.isLineElement(var, true), mLineElementAnnotator.isLineElement(var, false));
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
                writer.write("SampleId,ClusterId,ClusterCount,Id,Type,Ploidy");

                // position and copy number
                writer.write(",ChrStart,PosStart,OrientStart,ArmStart,AdjAFStart,AdjCNStart,AdjCNChgStart");

                writer.write(",ChrEnd,PosEnd,OrientEnd,ArmEnd,AdjAFEnd,AdjCNEnd,AdjCNChgEnd");

                // SV info
                writer.write(",Homology,InexactHOStart,InexactHOEnd,InsertSeq,Imprecise,PONCount,QualScore");

                // location attributes
                writer.write(",FSStart,FSEnd,LEStart,LEEnd,DupBEStart,DupBEEnd,ArmCountStart,ArmExpStart,ArmCountEnd,ArmExpEnd");

                // cluster-level info
                writer.write(",ClusterDesc,Consistency,ArmCount");

                // linked pair info
                writer.write(",LnkSvStart,LnkTypeStart,LnkLenStart,LnkSvEnd,LnkTypeEnd,LnkLenEnd");

                // GRIDDS caller info
                writer.write(",AsmbStart,AsmbEnd,AsmbMatchStart,AsmbMatchEnd");

                // chain info
                writer.write(",ChainId,ChainCount,ChainTICount,ChainDBCount,ChainIndex");

                // proximity info
                writer.write(",NearestLen,NearestType");

                // transitive info
                writer.write(",TransType,TransLen,TransSvLinks");

                // gene info
                // writer.write(",GeneStart,GeneDriverStart,GeneTypeStart,GeneEnd,GeneDriverEnd,GeneTypeEnd");

                writer.newLine();
            }

            int lineCount = 0;
            int svCount = 0;

            for(final SvCluster cluster : mClusters)
            {
                for (final SvClusterData var : cluster.getSVs())
                {
                    final StructuralVariantData dbData = var.getSvData();

                    ++svCount;

                    String typeStr = var.isNullBreakend() ? "SGL" : var.type().toString();
                    writer.write(
                            String.format("%s,%d,%d,%s,%s,%.2f",
                                    mSampleId, cluster.getId(), cluster.getCount(), var.id(), typeStr, dbData.ploidy()));

                    writer.write(
                            String.format(",%s,%d,%d,%s,%.2f,%.2f,%.2f,%s,%d,%d,%s,%.2f,%.2f,%.2f",
                                    var.chromosome(true), var.position(true), var.orientation(true), var.getStartArm(),
                                    dbData.adjustedStartAF(), dbData.adjustedStartCopyNumber(), dbData.adjustedStartCopyNumberChange(),
                                    var.chromosome(false), var.position(false), var.orientation(false), var.getEndArm(),
                                    dbData.adjustedEndAF(), dbData.adjustedEndCopyNumber(), dbData.adjustedEndCopyNumberChange()));

                    writer.write(
                            String.format(",%s,%d,%d,%s,%s,%d,%.0f",
                                    dbData.insertSequence().isEmpty() && var.type() != StructuralVariantType.INS ? dbData.homology() : "",
                                    dbData.inexactHomologyOffsetStart(), dbData.inexactHomologyOffsetEnd(),
                                    dbData.insertSequence(), dbData.imprecise(), var.getPonCount(), dbData.qualityScore()));

                    writer.write(
                            String.format(",%s,%s,%s,%s,%s,%s,%s",
                                    var.isStartFragileSite(), var.isEndFragileSite(),
                                    var.isStartLineElement(), var.isEndLineElement(),
                                    var.isDupBEStart(), var.isDupBEEnd(),
                                    mClusteringMethods.getChrArmData(var)));

                    writer.write(
                            String.format(",%s,%d,%d",
                                    cluster.getDesc(), cluster.getConsistencyCount(), cluster.getChromosomalArmCount()));

                    // linked pair info
                    final SvLinkedPair startLP = cluster.getLinkedPair(var, true);
                    String startLinkStr = "0,,-1";
                    String assemblyMatchStart = var.getAssemblyMatchType(true);
                    if(startLP != null)
                    {
                        startLinkStr = String.format("%s,%s,%d",
                                startLP.first().equals(var) ? startLP.second().id() : startLP.first().id(), startLP.linkType(), startLP.length());
                    }

                    final SvLinkedPair endLP = cluster.getLinkedPair(var, false);
                    String endLinkStr = "0,,-1";
                    String assemblyMatchEnd = var.getAssemblyMatchType(false);
                    if(endLP != null)
                    {
                        endLinkStr = String.format("%s,%s,%d",
                                endLP.first().equals(var) ? endLP.second().id() : endLP.first().id(), endLP.linkType(), endLP.length());
                    }

                    if(assemblyMatchStart.equals(ASSEMBLY_MATCH_ASMB_ONLY) || assemblyMatchEnd.equals(ASSEMBLY_MATCH_ASMB_ONLY))
                    {
                        LOGGER.debug("sample({}) var({}) has unlinked assembly TIs", mSampleId, var.posId());
                    }

                    // assembly info
                    writer.write(String.format(",%s,%s,%s,%s,%s,%s",
                            startLinkStr, endLinkStr, var.getAssemblyStart(), var.getAssemblyEnd(), assemblyMatchStart, assemblyMatchEnd));

                    // chain info
                    final SvChain chain = cluster.findChain(var);
                    String chainStr = ",0,0,0,0,0";

                    if(chain != null)
                    {
                        chainStr = String.format(",%d,%d,%d,%d,%d",
                                chain.getId(), chain.getLinkCount(), chain.getTICount(), chain.getDBCount(), chain.getSvIndex(var));
                    }

                    writer.write(chainStr);

                    writer.write(String.format(",%d,%s", var.getNearestSvDistance(), var.getNearestSvRelation()));

                    writer.write(String.format(",%s,%d,%s", var.getTransType(), var.getTransLength(), var.getTransSvLinks()));

                    /*
                    final SvGeneData geneStart = var.getStartGeneData();
                    if(geneStart != null)
                    {
                        writer.write(String.format(",%s,%s,%s", geneStart.gene(), geneStart.driver(), geneStart.driverType()));
                    }
                    else
                    {
                        writer.write(String.format(",,,"));
                    }

                    final SvGeneData geneEnd = var.getEndGeneData();
                    if(geneEnd != null)
                    {
                        writer.write(String.format(",%s,%s,%s", geneEnd.gene(), geneEnd.driver(), geneEnd.driverType()));
                    }
                    else
                    {
                        writer.write(String.format(",,,"));
                    }
                    */

                    ++lineCount;
                    writer.newLine();

                    if(svCount != lineCount)
                    {
                        LOGGER.error("inconsistent output");
                    }
                }
            }

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing to outputFile: {}", e.toString());
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
        mPc4.logStats(false);
        mPc5.logStats(false);
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

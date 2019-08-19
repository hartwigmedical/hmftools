package com.hartwig.hmftools.linx.drivers;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.UNKNOWN;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.findCentromereBreakendIndex;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CN_SEG_DATA_CN_BEFORE;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.GAIN;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.GAIN_ARM;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.GAIN_CHR;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_ARM;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_CHR;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_ARM_SV;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_CENTRO_SV;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_DEL;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CENTROMERE_CN;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.P_ARM_TELOMERE_CN;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.Q_ARM_TELOMERE_CN;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriverFile;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DriverGeneAnnotator
{
    private final DatabaseAccess mDbAccess;
    private SvGeneTranscriptCollection mGeneTransCache;
    private final LinxConfig mConfig;

    private List<DriverCatalog> mDriverCatalog;
    private List<GeneCopyNumber> mGeneCopyNumberData;
    private BufferedWriter mFileWriter;
    private String mOutputDir;
    private double mSamplePloidy;
    private List<DriverGeneData> mDriverGeneDataList;

    private List<LinxDriver> mDriverOutputList;

    private PerformanceCounter mPerfCounter;

    // references only
    private String mSampleId;
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private CnDataLoader mCopyNumberData;
    private Map<String, List<GeneCopyNumber>> mSampleGeneCopyNumberMap; // loaded from file to avoid DB hits on the massive table
    private VisualiserWriter mVisWriter;

    private static final String GCN_DATA_FILE = "gcn_data_file";

    private static final Logger LOGGER = LogManager.getLogger(DriverGeneAnnotator.class);

    public DriverGeneAnnotator(DatabaseAccess dbAccess, SvGeneTranscriptCollection geneTranscriptCollection,
            final LinxConfig config, CnDataLoader cnDataLoader)
    {
        mDbAccess = dbAccess;
        mGeneTransCache = geneTranscriptCollection;
        mConfig = config;
        mCopyNumberData = cnDataLoader;

        mDriverCatalog = Lists.newArrayList();
        mGeneCopyNumberData = Lists.newArrayList();
        mDriverGeneDataList = Lists.newArrayList();
        mDriverOutputList = Lists.newArrayList();
        mSampleGeneCopyNumberMap = new HashMap();

        mFileWriter = null;
        mOutputDir = mConfig.OutputDataPath;
        mSamplePloidy = 0;

        mGeneTransCache.createGeneNameIdMap();
        mVisWriter = null;

        mPerfCounter = new PerformanceCounter("Drivers");
    }

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(GCN_DATA_FILE, true, "Cache of driver-matched gene copy number data");
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        initialiseGeneData(cmd.getOptionValue(GCN_DATA_FILE,""));

        return true;
    }

    public void setSamplePloidy(double ploidy)
    {
        mSamplePloidy = ploidy;
    }

    public void setVisWriter(VisualiserWriter writer) { mVisWriter = writer; }
    public final List<DriverGeneData> getDriverGeneDataList() { return mDriverGeneDataList; }

    private void initialiseGeneData(final String geneCopyNumberFile)
    {
        loadGeneCopyNumberDataFile(geneCopyNumberFile);
    }

    private void loadDriverCatalog(final String sampleId)
    {
        mDriverCatalog.clear();
        mDriverGeneDataList.clear();
        mDriverCatalog.addAll(mDbAccess.readDriverCatalog(sampleId));

        LOGGER.debug("retrieved {} driver gene records", mDriverCatalog.size());
    }

    private void loadGeneCopyNumberData(final String sampleId)
    {
        mGeneCopyNumberData.clear();

        if(!mSampleGeneCopyNumberMap.isEmpty())
        {
            final List<GeneCopyNumber> gcnDataList = mSampleGeneCopyNumberMap.get(sampleId);

            if(gcnDataList != null)
                mGeneCopyNumberData.addAll(gcnDataList);

            return;
        }

        final List<String> driverGenes = mDriverCatalog.stream().map(x -> x.gene()).collect(Collectors.toList());

        mGeneCopyNumberData = mDbAccess.readGeneCopynumbers(sampleId, driverGenes);
    }

    public void annotateSVs(final String sampleId, final List<SvCluster> clusters, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        if(mDbAccess == null)
        {
            LOGGER.error("driver analysis requires DB connection");
            return;
        }

        mPerfCounter.start();

        final PurityContext purityContext = mDbAccess.readPurityContext(sampleId);

        if(purityContext != null)
            setSamplePloidy(purityContext.bestFit().ploidy());

        mSampleId = sampleId;
        mChrBreakendMap = chrBreakendMap;

        loadDriverCatalog(sampleId);

        if (mDriverCatalog.isEmpty())
            return;

        loadGeneCopyNumberData(sampleId);

        // types handled:
        // - TSG biallelic point mutations
        // - AMPs
        // - DELs

        for (final DriverCatalog driverGene : mDriverCatalog)
        {
            if(!driverTypeHandled(driverGene))
                continue;

            final EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(driverGene.gene());

            if (geneData == null)
            {
                LOGGER.warn("driver gene({}) data not found", driverGene.gene());
                continue;
            }

            final TranscriptData canonicalTrans = mGeneTransCache.getTranscriptData(geneData.GeneId, "");

            GeneCopyNumber geneCN = findGeneCopyNumber(driverGene);

            DriverGeneData dgData = new DriverGeneData(driverGene, geneData, canonicalTrans, geneCN);
            mDriverGeneDataList.add(dgData);

            if(dgData.DriverData.category() == TSG)
            {
                if(dgData.DriverData.driver() == DriverType.DEL)
                {
                    annotateDeleteEvent(dgData);
                }
                else
                {
                    annotateBiallelicEvent(dgData);
                }
            }
            else if (dgData.DriverData.category() == ONCO)
            {
                final List<SvBreakend> breakendList = mChrBreakendMap.get(dgData.GeneData.Chromosome);
                annotateAmplification(dgData, breakendList);
            }

            writeDriverData(dgData);
        }

        cacheSampleDriverData();

        mChrBreakendMap = null;

        mPerfCounter.stop();
    }

    private boolean driverTypeHandled(final DriverCatalog driverGene)
    {
        if(driverGene.category() == TSG && driverGene.driver() == DriverType.DEL)
            return true;

        if(driverGene.category() == TSG && driverGene.biallelic())
        {
            // need to look for an LOH
            return true;
        }

        if(driverGene.category() == ONCO && driverGene.driver() == DriverType.AMP)
            return true;

        return false;
    }

    private void annotateDeleteEvent(final DriverGeneData dgData)
    {
        if (dgData.GeneCN == null)
        {
            LOGGER.warn("sample({}) gene({}) min copy number data not found for DEL driver", mSampleId, dgData.GeneData.GeneName);
            return;
        }

        // find the LOH and hom-loss events which caused this DEL
        long minRegionStart = dgData.GeneCN.minRegionStart();
        long minRegionEnd = dgData.GeneCN.minRegionEnd();

        for(final LohEvent lohEvent : mCopyNumberData.getLohData())
        {
            if(!lohEvent.Chromosome.equals(dgData.GeneData.Chromosome))
                continue;

            // the LOH needs to straddle one or the other of the min-gene breakends
            if(lohEvent.PosStart > minRegionEnd || lohEvent.PosEnd < minRegionStart)
                continue;

            LOGGER.debug("gene({}) minCnRegion({} -> {}) covered by LOH({})",
                    dgData.GeneData.GeneName, minRegionStart, minRegionEnd, lohEvent);

            if((lohEvent.PosStart < minRegionEnd && lohEvent.PosEnd > minRegionStart) || !lohEvent.doubleSvEvent())
            {
                // LOH covers the whole gene or extends to end of arm (ie not just 2 SVs
                if(lohEvent.doubleSvEvent())
                {
                    DriverGeneEvent event = new DriverGeneEvent(LOH);
                    event.setLohEvent(lohEvent);

                    // one or both breakendscan be null if not matched
                    event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_DEL);
                    dgData.addEvent(event);
                }
                else if(lohEvent.armLoss())
                {
                    dgData.addEvent(new DriverGeneEvent(LOH_ARM));
                }
                else if(lohEvent.chromosomeLoss())
                {
                    dgData.addEvent(new DriverGeneEvent(LOH_CHR));
                }
                else if(lohEvent.isSvEvent())
                {
                    // call SV + rest of arm loss an LOH as well
                    DriverGeneEvent event = new DriverGeneEvent(LOH);
                    event.setLohEvent(lohEvent);
                    event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_ARM_SV);
                    dgData.addEvent(event);
                }
            }
            else
            {
                // cannot assign one event as the LOH and one as the DEL, so assign both as DELs
                DriverGeneEvent event = new DriverGeneEvent(DriverEventType.DEL);
                event.setLohEvent(lohEvent);

                // one or both breakendscan be null if not matched
                event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_DEL);
                dgData.addEvent(event);
            }

            for(final HomLossEvent homLoss : lohEvent.getHomLossEvents())
            {
                // allow there to be more than one?
                if(homLoss.PosStart <= minRegionEnd && homLoss.PosEnd >= minRegionStart)
                {
                    LOGGER.debug("gene({}) minCnRegion({} -> {}) covered by hom-loss({})",
                            dgData.GeneData.GeneName, minRegionStart, minRegionEnd, homLoss);

                    // DEL covers the whole gene or extends to end of arm
                    DriverGeneEvent event = new DriverGeneEvent(DriverEventType.DEL);
                    event.setHomLossEvent(homLoss);

                    // one or both breakends can be null if not matched
                    event.addSvBreakendPair(homLoss.getBreakend(true), homLoss.getBreakend(false), SV_DRIVER_TYPE_DEL);
                    dgData.addEvent(event);
                    break;
                }
            }

            break;
        }
    }

    private void annotateBiallelicEvent(DriverGeneData dgData)
    {
        // look for an LOH covering any part of the coding region
        if(dgData.TransData.CodingStart == null || dgData.TransData.CodingEnd == null)
            return;

        long codingStart = dgData.TransData.CodingStart;
        long codingEnd = dgData.TransData.CodingEnd;

        for(final LohEvent lohEvent : mCopyNumberData.getLohData())
        {
            if(!lohEvent.Chromosome.equals(dgData.GeneData.Chromosome))
                continue;

            if(lohEvent.PosStart > codingEnd || lohEvent.PosEnd < codingStart)
                continue;

            LOGGER.debug("gene({}) coding region({} -> {}) covered by LOH({})",
                    dgData.GeneData.GeneName, codingStart, codingEnd, lohEvent);

            if(lohEvent.doubleSvEvent())
            {
                DriverGeneEvent event = new DriverGeneEvent(LOH);
                event.setLohEvent(lohEvent);

                // one or both breakendscan be null if not matched
                event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_DEL);
                dgData.addEvent(event);
            }
            else if(lohEvent.armLoss())
            {
                dgData.addEvent(new DriverGeneEvent(LOH_ARM));
            }
            else if(lohEvent.chromosomeLoss())
            {
                dgData.addEvent(new DriverGeneEvent(LOH_CHR));
            }
            else if(lohEvent.isSvEvent())
            {
                DriverGeneEvent event = new DriverGeneEvent(LOH_ARM);
                event.setLohEvent(lohEvent);
                event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_ARM_SV);
                dgData.addEvent(event);
            }

            break;
        }
    }

    private void annotateAmplification(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        if(mSamplePloidy == 0)
            return;

        LOGGER.debug("gene({}) chromosome({}) position({} -> {}) minCN({})",
                dgData.GeneData.GeneName, dgData.GeneData.Chromosome, dgData.TransData.TransStart, dgData.TransData.TransEnd,
                formatPloidy(dgData.GeneCN.minCopyNumber()));

        checkArmAmplification(dgData);
        checkClusterAmplification(dgData, breakendList, true);

        if(dgData.getEvents().isEmpty())
        {
            checkClusterAmplification(dgData, breakendList, false);
        }

        if(dgData.getEvents().isEmpty())
        {
            checkAmpInCentromere(dgData, breakendList);
        }

        if(dgData.getEvents().isEmpty())
        {
            LOGGER.debug("gene({}) AMP gain no event found", dgData.GeneData.GeneName);

            // otherwise no event
            DriverGeneEvent event = new DriverGeneEvent(GAIN);
            dgData.addEvent(event);
        }
    }

    private void checkArmAmplification(final DriverGeneData dgData)
    {
        // check for an arm or whole chromosome amplified above the sample ploidy
        double[] cnData = mCopyNumberData.getChrCopyNumberMap().get(dgData.GeneData.Chromosome);

        // find the copy number leading into the centromere on the arm for this gene
        double[] centromereCnData = mCopyNumberData.getCentromereCopyNumberData(
                dgData.GeneData.Chromosome, dgData.Arm == CHROMOSOME_ARM_P);

        double centromereCopyNumber = centromereCnData[CN_SEG_DATA_CN_BEFORE];

        double telomereCopyNumber = min(cnData[P_ARM_TELOMERE_CN], cnData[Q_ARM_TELOMERE_CN]);

        double armCopyNumber = min(centromereCopyNumber,
                dgData.Arm == CHROMOSOME_ARM_P ? cnData[P_ARM_TELOMERE_CN] : cnData[Q_ARM_TELOMERE_CN]);

        double chromosomeCopyNumber = min(centromereCopyNumber, telomereCopyNumber);

        if(chromosomeCopyNumber/mSamplePloidy > 2)
        {
            LOGGER.debug("gene({}) AMP gain from chromosomeCN({}) vs samplePloidy({})",
                    dgData.GeneData.GeneName, formatPloidy(chromosomeCopyNumber), formatPloidy(mSamplePloidy));

            dgData.addEvent(new DriverGeneEvent(GAIN_CHR));
        }
        else if(armCopyNumber/mSamplePloidy > 2)
        {
            LOGGER.debug("gene({}) AMP gain from armCN({}) vs samplePloidy({})",
                    dgData.GeneData.GeneName, formatPloidy(armCopyNumber), formatPloidy(mSamplePloidy));

            dgData.addEvent(new DriverGeneEvent(GAIN_ARM));
        }
    }

    private DriverAmpData checkClusterForAmplification(
            final DriverGeneData dgData, final List<SvBreakend> breakendList, final SvBreakend startBreakend, boolean traverseUp,
            Map<SvBreakend,Double> opposingBreakends)
    {
        long transStart = dgData.TransData.TransStart;
        long transEnd = dgData.TransData.TransEnd;

        double startCopyNumber = startBreakend.getCopyNumber(traverseUp);

        final SvCluster targetCluster = startBreakend.getCluster();
        final List<SvBreakend> clusterBreakends = targetCluster.getChrBreakendMap().get(dgData.GeneData.Chromosome);
        int clusterStartIndex = clusterBreakends.get(0).getChrPosIndex();
        int clusterEndIndex = clusterBreakends.get(clusterBreakends.size() - 1).getChrPosIndex();

        boolean inSegment = true;
        double segmentStartCopyNumber = startCopyNumber;
        double netClusterCNChange = 0;
        int segmentCount = 0;
        int breakendCount = 0;

        int index = startBreakend.getChrPosIndex();
        SvBreakend breakend = null;
        SvBreakend segStartBreakend = startBreakend;

        while (true)
        {
            if (breakend != null)
                index += traverseUp ? 1 : -1;

            if(index < 0 || index > breakendList.size())
                break;

            breakend = breakendList.get(index);

            if ((traverseUp && breakend.position() > transStart) || (!traverseUp && breakend.position() < transEnd))
                break;

            final SvCluster cluster = breakend.getCluster();

            if(cluster != targetCluster)
            {
                if(inSegment)
                {
                    // record details of this segment
                    double endCopyNumber = breakend.getCopyNumber(traverseUp);
                    double segmentCNChange = endCopyNumber - segmentStartCopyNumber;
                    netClusterCNChange += segmentCNChange;

                    LOGGER.trace("gene({}) cluster({}) adding segment CN({} -> {} chg={}) net({}) breakends({} -> {})",
                            dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(segmentStartCopyNumber), formatPloidy(endCopyNumber),
                            formatPloidy(segmentCNChange), formatPloidy(netClusterCNChange), breakend, segStartBreakend);

                    // if this next breakend lowers

                    inSegment = false;
                }

                boolean isOpposingBreakend = (breakend.orientation() == 1) == traverseUp;

                if(netClusterCNChange > 0 && isOpposingBreakend)
                {
                    Double remainingCNChange = opposingBreakends.get(breakend);

                    if(remainingCNChange == null)
                    {
                        remainingCNChange = breakend.copyNumberChange();
                    }

                    if(remainingCNChange > netClusterCNChange)
                    {
                        LOGGER.trace("gene({}) cluster({}) netCN({}) cancelled by breakend({}) cnChange(orig={} remain={})",
                                dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(netClusterCNChange),
                                breakend, formatPloidy(breakend.copyNumberChange()), formatPloidy(remainingCNChange));

                        remainingCNChange -= netClusterCNChange;
                        netClusterCNChange = 0;

                        opposingBreakends.put(breakend, remainingCNChange);
                    }
                    else
                    {
                        LOGGER.trace("gene({}) cluster({}) netCN({}) reducing by breakend({}) cnChange({})",
                                dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(netClusterCNChange),
                                breakend, formatPloidy(breakend.copyNumberChange()));

                        netClusterCNChange -= remainingCNChange;
                        opposingBreakends.put(breakend, 0.0); // keep so it won't be registered again but zero out
                    }
                }

                // break if there are no more breakends in this cluster
                // if(index < clusterStartIndex || index > clusterEndIndex)
                //    break;
            }
            else if(cluster == targetCluster)
            {
                ++breakendCount;

                if(!inSegment)
                {
                    segmentStartCopyNumber = breakend.getCopyNumber(traverseUp);
                    segStartBreakend = breakend;
                    inSegment = true;
                }
            }
        }

        double geneMinCopyNumber = dgData.GeneCN.minCopyNumber();

        double clusterCNChange = netClusterCNChange;

        if(inSegment)
        {
            clusterCNChange += geneMinCopyNumber - segmentStartCopyNumber;

            LOGGER.trace("gene({}) cluster({}) open segment startCN({}) net({}) start breakend({})",
                    dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(segmentStartCopyNumber), formatPloidy(clusterCNChange),
                    segStartBreakend);
        }

        if(clusterCNChange > 0 && !copyNumbersEqual(clusterCNChange - startCopyNumber, startCopyNumber))
        {
            LOGGER.debug("gene({}) cluster({}) copy number gain({}) vs startCN({}) segments({}: CN={}) traversal({})",
                    dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(clusterCNChange), formatPloidy(startCopyNumber),
                    segmentCount, formatPloidy(netClusterCNChange), traverseUp ? "up" : "down");

            return new DriverAmpData(targetCluster, traverseUp, breakendCount, segmentCount, startCopyNumber, clusterCNChange);
        }

        return null;
    }

    private void checkClusterAmplification(final DriverGeneData dgData, final List<SvBreakend> breakendList, boolean restrictToArm)
    {
        if (breakendList == null || breakendList.isEmpty())
            return;

        // sum up breakend ploidies from telomere to centromere for the gene in question net off ploidy within a cluster
        long transStart = dgData.TransData.TransStart;
        long transEnd = dgData.TransData.TransEnd;

        int centromereIndex = findCentromereBreakendIndex(breakendList, dgData.Arm);

        int startIndex = (dgData.Arm == CHROMOSOME_ARM_P || !restrictToArm) ? 0 : centromereIndex;
        int endIndex = (dgData.Arm == CHROMOSOME_ARM_P && restrictToArm ) ? centromereIndex : breakendList.size() - 1;

        Map<SvCluster,DriverAmpData> clusterAmpData = Maps.newHashMap();

        for(int i = 0; i <= 1; ++i)
        {
            boolean traverseUp = (i == 0);

            List<SvCluster> processedClusters = Lists.newArrayList();
            Map<SvBreakend,Double> opposingBreakends = Maps.newHashMap();

            int index = traverseUp ? startIndex : endIndex;
            SvBreakend breakend = null;

            while (true)
            {
                if (breakend != null)
                    index += traverseUp ? 1 : -1;

                if (index < startIndex || index > endIndex)
                    break;

                breakend = breakendList.get(index);

                // make note of any breakend preceding the gene in the direction of traversal
                if ((traverseUp && breakend.position() > transEnd) || (!traverseUp && breakend.position() < transStart))
                    break;

                final SvCluster cluster = breakend.getCluster();

                if (cluster.hasLinkingLineElements() || processedClusters.contains(cluster))
                    continue;

                processedClusters.add(cluster);

                // proceeed from this point until the start of the gene
                DriverAmpData ampData = checkClusterForAmplification(dgData, breakendList, breakend, traverseUp, opposingBreakends);

                if(ampData == null)
                    continue;

                DriverAmpData existingAmpData = clusterAmpData.get(cluster);

                if(existingAmpData == null || existingAmpData.NetCNChange < ampData.NetCNChange)
                {
                    clusterAmpData.put(cluster, ampData);
                }
            }
        }

        if(clusterAmpData.isEmpty())
            return;

        // from the identified AMP clusters, find the max and percentage contributions of each
        double maxCnChange = clusterAmpData.values().stream().mapToDouble(x -> x.NetCNChange).max().getAsDouble();

        for(Map.Entry<SvCluster,DriverAmpData> entry : clusterAmpData.entrySet())
        {
            final SvCluster cluster = entry.getKey();
            final DriverAmpData ampData = entry.getValue();

            if(ampData.NetCNChange / maxCnChange < 0.1)
                continue;

            LOGGER.debug("gene({}) cluster({}) adding AMP data: {}",
                    dgData.GeneData.GeneName, cluster.id(), ampData);

            DriverGeneEvent event = new DriverGeneEvent(GAIN);
            event.setCluster(cluster);
            event.setAmpData(ampData);

            event.setSvInfo(String.format("%s;%d;%d;%.1f;%.1f",
                    ampData.TraverseUp, ampData.BreakendCount, ampData.SegmentCount, ampData.StartCopyNumber, ampData.NetCNChange));
            dgData.addEvent(event);
        }
    }

    private void checkAmpInCentromere(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        // check for a CN gain across the centromere
        int centromereIndex = findCentromereBreakendIndex(breakendList, dgData.Arm);
        int preCentromereIndex = dgData.Arm == CHROMOSOME_ARM_P ?  centromereIndex : centromereIndex - 1;
        int postCentromereIndex = dgData.Arm == CHROMOSOME_ARM_P ? centromereIndex + 1 : centromereIndex;

        if(preCentromereIndex >= 0 && postCentromereIndex < breakendList.size())
        {
            final SvBreakend preCentroBreakend = breakendList.get(preCentromereIndex);
            final SvBreakend postCentroBreakend = breakendList.get(postCentromereIndex);

            double preCopyNumber = preCentroBreakend.getCopyNumber(false);
            double postCopyNumber = postCentroBreakend.getCopyNumber(true);
            double centromereCNChange =
                    dgData.Arm == CHROMOSOME_ARM_P ? preCopyNumber - postCopyNumber : postCopyNumber - preCopyNumber;

            if (centromereCNChange > 0 && centromereCNChange > mSamplePloidy * 2)
            {
                LOGGER.debug("gene({}) copy number change across centromere({} -> {} = {})",
                        dgData.GeneData.GeneName, formatPloidy(preCopyNumber), formatPloidy(postCopyNumber),
                        formatPloidy(centromereCNChange));

                DriverGeneEvent event = new DriverGeneEvent(GAIN);
                event.setSvInfo(SV_DRIVER_TYPE_CENTRO_SV);
                dgData.addEvent(event);
            }
        }
    }

    private final GeneCopyNumber findGeneCopyNumber(final DriverCatalog driverGene)
    {
        for(GeneCopyNumber geneCN : mGeneCopyNumberData)
        {
            if(driverGene.gene().equals(geneCN.gene()))
                return geneCN;
        }

        return null;
    }

    private void cacheSampleDriverData()
    {
        if(mConfig.hasMultipleSamples())
            return;

        if(mDriverOutputList.isEmpty())
            return;

        if(mConfig.UploadToDB)
        {
            mDbAccess.writeSvDrivers(mSampleId, mDriverOutputList);
        }

        try
        {
            final String driversFile = LinxDriverFile.generateFilename(mOutputDir, mSampleId);
            LinxDriverFile.write(driversFile, mDriverOutputList);
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write drivers file: {}", e.toString());
        }
    }

    private void writeDriverData(final DriverGeneData dgData)
    {
        // convert to a sample driver record
        for(final DriverGeneEvent driverEvent : dgData.getEvents())
        {
            int clusterId = driverEvent.getCluster() != null ? driverEvent.getCluster().id() : -1;

            mDriverOutputList.add(ImmutableLinxDriver.builder()
                    .clusterId(clusterId)
                    .gene(dgData.GeneData.GeneName)
                    .eventType(driverEvent.Type.toString())
                    .build());

            mVisWriter.addGeneExonData(clusterId, dgData.GeneData.GeneId, dgData.GeneData.GeneName,
                    "", 0, dgData.GeneData.Chromosome, "DRIVER");
        }

        if(!mConfig.hasMultipleSamples())
            return;

        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += "SVA_DRIVERS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Gene,Category,DriverType,LikelihoodMethod");
                mFileWriter.write(",FullyMatched,EventType,ClusterId,ClusterCount,ResolvedType");
                mFileWriter.write(",PosStart,PosEnd,SvIdStart,SvIdEnd,SvMatchType,SvPloidy");
                mFileWriter.write(",SamplePloidy,Chromosome,Arm,GeneMinCN,CentromereCN,TelomereCN");
                mFileWriter.newLine();
            }

            BufferedWriter writer = mFileWriter;

            final DriverCatalog driverGene = dgData.DriverData;
            final EnsemblGeneData geneData = dgData.GeneData;
            double[] cnData = mCopyNumberData.getChrCopyNumberMap().get(geneData.Chromosome);

            for(final DriverGeneEvent driverEvent : dgData.getEvents())
            {
                final SvBreakend[] breakendPair = driverEvent.getBreakendPair();

                writer.write(String.format("%s,%s,%s,%s,%s,%s,%s",
                        mSampleId, driverGene.gene(), driverGene.category(), driverGene.driver(), driverGene.likelihoodMethod(),
                        dgData.fullyMatched(), driverEvent.Type));

                final SvCluster cluster = driverEvent.getCluster();

                if(cluster != null)
                {
                    writer.write(String.format(",%d,%d,%s", cluster.id(), cluster.getSvCount(), cluster.getResolvedType()));
                }
                else
                {
                    writer.write(String.format(",-1,0,"));
                }

                // breakend info if present

                long posStart = breakendPair[SE_START] != null ? breakendPair[SE_START].position() : 0;
                long posEnd = breakendPair[SE_END] != null ? breakendPair[SE_END].position() : 0;
                String svIdStart = breakendPair[SE_START] != null ? breakendPair[SE_START].getSV().idStr() : "-1";
                String svIdEnd = breakendPair[SE_END] != null ? breakendPair[SE_END].getSV().idStr() : "-1";

                double svPloidy = breakendPair[SE_START] != null ? breakendPair[SE_START].ploidy()
                        : (breakendPair[SE_END] != null ? breakendPair[SE_END].ploidy() : 0);

                writer.write(String.format(",%d,%d,%s,%s,%s,%.2f",
                        posStart, posEnd, svIdStart, svIdEnd, driverEvent.getSvInfo(), svPloidy));

                writer.write(String.format(",%.2f,%s,%s,%.2f,%.2f,%.2f",
                        mSamplePloidy, geneData.Chromosome, dgData.Arm, dgData.GeneCN != null ? dgData.GeneCN.minCopyNumber() : -1,
                        cnData[CENTROMERE_CN], dgData.Arm == CHROMOSOME_ARM_P ? cnData[P_ARM_TELOMERE_CN] : cnData[Q_ARM_TELOMERE_CN]));

                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing cluster-data to outputFile: {}", e.toString());
        }
    }

    public void close()
    {
        mPerfCounter.logStats();

        closeBufferedWriter(mFileWriter);
    }

    private static int GENE_CN_DATA_FILE_ITEM_COUNT = 6;

    private void loadGeneCopyNumberDataFile(final String gcnFileName)
    {
        if(gcnFileName.isEmpty() || !Files.exists(Paths.get(gcnFileName)))
            return;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(gcnFileName));

            // skip field names
            String line = fileReader.readLine();

            if (line == null)
            {
                LOGGER.error("Empty gene copy number CSV file({})", gcnFileName);
                return;
            }

            String currentSample = "";
            List<GeneCopyNumber> gcnDataList = null;
            int rowCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if (items.length != GENE_CN_DATA_FILE_ITEM_COUNT)
                {
                    LOGGER.error("GCN file invalid item count({}) vs expected({})", items.length, GENE_CN_DATA_FILE_ITEM_COUNT);
                    break;
                }

                String sampleId = items[0];

                if(currentSample.isEmpty() || !currentSample.equals(sampleId))
                {
                    gcnDataList = Lists.newArrayList();
                    currentSample = sampleId;
                    mSampleGeneCopyNumberMap.put(currentSample, gcnDataList);
                }

                // sampleId,minCopyNumber,minRegionStart,minRegionEnd,minMinorAllelePloidy

                int index = 1;

                GeneCopyNumber gcnData = ImmutableGeneCopyNumber.builder()
                        .gene(items[index++])
                        .minCopyNumber(Double.parseDouble(items[index++]))
                        .minRegionStart(Integer.parseInt(items[index++]))
                        .minRegionEnd(Integer.parseInt(items[index++]))
                        .minMinorAllelePloidy(Double.parseDouble(items[index++]))
                        .maxCopyNumber(0)
                        .somaticRegions(0)
                        .germlineHet2HomRegions(0)
                        .germlineHomRegions(0)
                        .minRegions(0)
                        .minRegionStartSupport(UNKNOWN)
                        .minRegionEndSupport(UNKNOWN)
                        .minRegionMethod(CopyNumberMethod.UNKNOWN)
                        .transcriptID("")
                        .transcriptVersion(0)
                        .chromosomeBand("")
                        .chromosome("")
                        .start(0)
                        .end(0)
                        .build();

                ++rowCount;

                gcnDataList.add(gcnData);
            }

            LOGGER.info("loaded {} gene copy-number records from file", rowCount);
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read gene copy number CSV file({}): {}", gcnFileName, e.toString());
        }
    }

    /* OLD AMP SV identification

            /*
        for(Map.Entry<SvCluster,List<SvBreakend>> entry : clusterBreakendsMap.entrySet())
        {
            final List<SvBreakend> clusterBreakends = entry.getValue();

            List<SvLinkedPair> spanningLinks = Lists.newArrayList();
            List<SvBreakend> facingFoldbacks = Lists.newArrayList();

            SvBreakend maxFacingPloidyBreakend = null;
            double maxFacingPloidy = 0;

            for(final SvBreakend breakend : clusterBreakends)
            {
                boolean breakendFacesGene = (breakend.position() < transStart && breakend.orientation() == -1)
                        || (breakend.position() > transEnd && breakend.orientation() == 1);

                if (!breakendFacesGene)
                    continue;

                if(breakend.isFoldback())
                {
                    boolean hasFoldback = facingFoldbacks.stream().anyMatch(x -> x.getSV() == breakend.getSV()
                            || x == breakend || x == breakend.getFoldbackBreakend());

                    if(!hasFoldback)
                        facingFoldbacks.add(breakend);
                }

                // make note of any TI which span the gene
                for(final SvLinkedPair pair : breakend.getSV().getLinkedPairs(breakend.usesStart()))
                {
                    if(spanningLinks.contains(pair))
                        continue;

                    if(pair.getBreakend(true).position() < transStart && pair.getBreakend(false).position() > transEnd)
                        spanningLinks.add(pair);
                }

                double breakendPloidy = max(breakend.ploidy(), breakend.copyNumberChange());

                if(breakendPloidy > 1 && (maxFacingPloidyBreakend == null || breakendPloidy > maxFacingPloidyBreakend.ploidy()))
                {
                    maxFacingPloidyBreakend = breakend;
                    maxFacingPloidy = breakendPloidy;
                }
            }

            // work out net and max ploidy in each direction
            double maxNetPloidy = 0;

            for(int i = 0; i <= 1; ++i)
            {
                boolean traverseUp = (i == 0);

                double netPloidy = 0;

                int index = traverseUp ? 0 : clusterBreakends.size() - 1;
                SvBreakend breakend = null;

                while (true)
                {
                    if (breakend != null)
                        index += traverseUp ? 1 : -1;

                    if (index < 0 || index >= clusterBreakends.size())
                        break;

                    breakend = breakendList.get(index);

                    boolean breakendFacesGene = ((breakend.orientation() != 1) == traverseUp);

                    if((traverseUp && breakend.position() > transEnd) || (!traverseUp && breakend.position() < transStart))
                        break;

                    netPloidy += breakendFacesGene ? breakend.ploidy() : -breakend.ploidy();
                }

                maxNetPloidy = max(maxNetPloidy, netPloidy);
            }

            if(spanningLinks.isEmpty() && maxFacingPloidy < mSamplePloidy && maxNetPloidy < mSamplePloidy && facingFoldbacks.isEmpty())
                continue;

            final SvCluster cluster = entry.getKey();

            // report on type of SV which caused the amp
            if(cluster.getSvCount() == 1 && cluster.getSV(0).type() == DUP && maxFacingPloidy > mSamplePloidy)
            {
                final SvVarData var = cluster.getSV(0);
                final String matchType = cluster.hasAnnotation(CLUSTER_ANNOT_DM) ? SV_DRIVER_TYPE_DM : SV_DRIVER_TYPE_DUP;
                DriverGeneEvent event = new DriverGeneEvent(GAIN);
                event.setCluster(cluster);
                event.addSvBreakendPair(var.getBreakend(true), var.getBreakend(false), matchType);
                dgData.addEvent(event);

                LOGGER.debug("gene({}) cluster({}) net ploidy({}) from DUP({}) type({}) sameArm({})",
                        dgData.GeneData.GeneName, cluster.id(), formatPloidy(maxFacingPloidy), var.toString(), matchType, restrictToArm);

                continue;
            }

            double spanningLinkPloidy = 0;
            double maxLinkPloidy = 0;

            if(!spanningLinks.isEmpty())
            {
                // calculate the total ploidy for this link from all chains it is in
                for (final SvLinkedPair pair : spanningLinks)
                {
                    for (final SvChain chain : cluster.getChains())
                    {
                        int linkRepeats = (int) chain.getLinkedPairs().stream().filter(x -> x.matches(pair)).count();

                        if (linkRepeats > 0)
                        {
                            double linkPloidy = linkRepeats * chain.ploidy();
                            ;
                            spanningLinkPloidy += linkPloidy;
                            maxLinkPloidy = max(maxLinkPloidy, linkPloidy);
                        }
                    }
                }
            }

            DriverGeneEvent event = new DriverGeneEvent(GAIN);
            event.setCluster(cluster);

            // record results as Type;Count;MaxPloidy;TotalPloidy

            if(spanningLinkPloidy > mSamplePloidy)
            {
                String linkInfo = ampSvInfo(
                        SV_DRIVER_TYPE_TI, clusterBreakends.size(), spanningLinks.size(), maxLinkPloidy, spanningLinkPloidy, !restrictToArm);

                event.setSvInfo(linkInfo);

                LOGGER.debug("gene({}) cluster({}) spanning TIs({}) maxPloidy({}) totalPloidy({}) sameArm({})",
                        dgData.GeneData.GeneName, cluster.id(), spanningLinks.size(), formatPloidy(maxLinkPloidy),
                        formatPloidy(spanningLinkPloidy), restrictToArm);
                dgData.addEvent(event);
                continue;
            }

            if(!facingFoldbacks.isEmpty())
            {
                double maxFoldbackPloidy = facingFoldbacks.stream().mapToDouble(x -> x.ploidy()).max().getAsDouble();

                double foldbackTotalPloidy = max(
                        facingFoldbacks.stream().filter(x -> x.orientation() == 1).mapToDouble(x -> x.ploidy()).max().orElse(0),
                        facingFoldbacks.stream().filter(x -> x.orientation() == -1).mapToDouble(x -> x.ploidy()).max().orElse(0));

                if(maxFoldbackPloidy > mSamplePloidy || foldbackTotalPloidy > mSamplePloidy)
                {
                    LOGGER.debug("gene({}) cluster({}) foldbacks({}) maxPloidy({}) totalFacing({}) sameArm({})",
                            dgData.GeneData.GeneName, cluster.id(), facingFoldbacks.size(), formatPloidy(maxFoldbackPloidy),
                            formatPloidy(foldbackTotalPloidy), restrictToArm);

                    String linkInfo = ampSvInfo(SV_DRIVER_TYPE_FOLDBACK, clusterBreakends.size(), facingFoldbacks.size(),
                            maxFoldbackPloidy, foldbackTotalPloidy, !restrictToArm);

                    event.setSvInfo(linkInfo);
                    dgData.addEvent(event);
                    continue;
                }
            }

            if(maxNetPloidy > mSamplePloidy)
            {
                LOGGER.debug("gene({}) cluster({}) max net ploidy({}) sameArm({})",
                        dgData.GeneData.GeneName, cluster.id(), formatPloidy(maxNetPloidy), restrictToArm);

                String linkInfo = ampSvInfo(
                        SV_DRIVER_TYPE_NET_PLOIDY, clusterBreakends.size(), 0, maxFacingPloidy, maxNetPloidy, !restrictToArm);
                event.setSvInfo(linkInfo);
                dgData.addEvent(event);
                continue;
            }

            if(maxFacingPloidy > mSamplePloidy)
            {
                LOGGER.debug("gene({}) cluster({}) max facing breakend({}) ploidy({}) sameArm({})",
                        dgData.GeneData.GeneName, cluster.id(), maxFacingPloidyBreakend, formatPloidy(maxFacingPloidy), restrictToArm);

                String linkInfo = ampSvInfo(
                        SV_DRIVER_TYPE_MAX_PLOIDY, clusterBreakends.size(), 1, maxFacingPloidy, maxNetPloidy, !restrictToArm);
                event.setSvInfo(linkInfo);
                dgData.addEvent(event);
                continue;
            }

            LOGGER.debug("gene({}) cluster({}) has breakend({}) without meeting amp criteria",
                    dgData.GeneData.GeneName, cluster.id(), clusterBreakends.size());
        }
     */

}

package com.hartwig.hmftools.linx.drivers;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod.DEL;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.UNKNOWN;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.MAX_COPY_NUM_DIFF;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.MAX_COPY_NUM_DIFF_PERC;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.copyNumbersEqual;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.isShortArmChromosome;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.TOTAL_CN_LOSS;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.GAIN;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.GAIN_ARM;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.GAIN_CHR;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.HOM_DEL_DISRUPTION;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.HOM_DUP_DISRUPTION;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_ARM;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_CHR;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_SV_CENTRO;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_SV_TELO;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_ARM_SV;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_DEL;
import static com.hartwig.hmftools.linx.types.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter.GENE_TYPE_DRIVER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;
import com.hartwig.hmftools.common.fusion.Transcript;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriverFile;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.cn.TelomereCentromereCnData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DriverGeneAnnotator
{
    private final DatabaseAccess mDbAccess;
    private final EnsemblDataCache mGeneTransCache;
    private final LinxConfig mConfig;

    private List<DriverCatalog> mDriverCatalog;
    private List<GeneCopyNumber> mGeneCopyNumberData;
    private BufferedWriter mFileWriter;
    private String mOutputDir;
    private double mSamplePloidy;
    private List<DriverGeneData> mDriverGeneDataList;
    final List<String> mReportableDelGeneIds;

    private List<LinxDriver> mDriverOutputList;

    private PerformanceCounter mPerfCounter;

    // references only
    private String mSampleId;
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private CnDataLoader mCopyNumberData;
    private Map<String, List<GeneCopyNumber>> mSampleGeneCopyNumberMap; // loaded from file to avoid DB hits on the massive table
    private VisualiserWriter mVisWriter;

    private static final String GCN_DATA_FILE = "gcn_data_file";

    public DriverGeneAnnotator(DatabaseAccess dbAccess, EnsemblDataCache geneTranscriptCollection,
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

        mReportableDelGeneIds = Lists.newArrayList();
        Set<String> reportableDelGenes = CNADrivers.reportableGeneDeletions();

        for(String geneName : reportableDelGenes)
        {
            final EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(geneName);

            if(geneData != null)
                mReportableDelGeneIds.add(geneData.GeneId);
        }

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

    private void loadDataFromDatabase()
    {
        final PurityContext purityContext = mDbAccess.readPurityContext(mSampleId);

        if(purityContext != null)
            setSamplePloidy(purityContext.bestFit().ploidy());

        mDriverCatalog.clear();

        // add records but filter out any previously added by Linx
        mDriverCatalog.addAll(
                mDbAccess.readDriverCatalog(mSampleId).stream()
                .filter(x -> x.driver() != DriverType.HOM_DISRUPTION)
                .collect(Collectors.toList()));

        LNX_LOGGER.debug("retrieved {} driver gene records", mDriverCatalog.size());

        if (!mDriverCatalog.isEmpty())
        {
            mGeneCopyNumberData.clear();

            if(!mSampleGeneCopyNumberMap.isEmpty())
            {
                final List<GeneCopyNumber> gcnDataList = mSampleGeneCopyNumberMap.get(mSampleId);

                if(gcnDataList != null)
                    mGeneCopyNumberData.addAll(gcnDataList);

                return;
            }

            final List<String> driverGenes = mDriverCatalog.stream().map(x -> x.gene()).collect(Collectors.toList());

            mGeneCopyNumberData = mDbAccess.readGeneCopynumbers(mSampleId, driverGenes);
        }
    }

    private void loadDataFromFile()
    {
        try
        {
            final PurityContext purityContext = FittedPurityFile.read(mConfig.PurpleDataPath, mSampleId);
            setSamplePloidy(purityContext.bestFit().ploidy());

            mDriverCatalog.addAll(DriverCatalogFile.read(DriverCatalogFile.generateFilename(mConfig.PurpleDataPath, mSampleId)));

            mGeneCopyNumberData.addAll(
                    GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilenameForReading(mConfig.PurpleDataPath, mSampleId)));
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load driver catalog or purity context: {}", e.toString());
            return;
        }
    }

    public void annotateSVs(final String sampleId, final Map<String, List<SvBreakend>> chrBreakendMap)
    {
        mPerfCounter.start();

        mSampleId = sampleId;
        mChrBreakendMap = chrBreakendMap;

        mDriverGeneDataList.clear();

        if(mDbAccess != null)
        {
            loadDataFromDatabase();
        }
        else if(!mConfig.PurpleDataPath.isEmpty())
        {
            loadDataFromFile();
        }

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
                LNX_LOGGER.warn("driver gene({}) data not found", driverGene.gene());
                continue;
            }

            final TranscriptData canonicalTrans = mGeneTransCache.getTranscriptData(geneData.GeneId, "");

            GeneCopyNumber geneCN = findGeneCopyNumber(driverGene.gene());

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

        findDisruptiveDelDrivers();

        writeSampleDriverData();

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
            LNX_LOGGER.warn("sample({}) gene({}) min copy number data not found for DEL driver", mSampleId, dgData.GeneData.GeneName);
            return;
        }

        // find the LOH and hom-loss events which caused this DEL
        int minRegionStart = (int)dgData.GeneCN.minRegionStart();
        int minRegionEnd = (int)dgData.GeneCN.minRegionEnd();

        for(final LohEvent lohEvent : mCopyNumberData.getLohData())
        {
            if(!lohEvent.Chromosome.equals(dgData.GeneData.Chromosome))
                continue;

            // the LOH needs to straddle one or the other of the min-gene breakends
            if(lohEvent.PosStart > minRegionEnd || lohEvent.PosEnd < minRegionStart)
                continue;

            LNX_LOGGER.debug("gene({}) minCnRegion({} -> {}) covered by LOH({})",
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
                    DriverEventType eventType = lohEvent.telomereLoss() ? LOH_SV_TELO : LOH_SV_CENTRO;
                    DriverGeneEvent event = new DriverGeneEvent(eventType);
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
                    LNX_LOGGER.debug("gene({}) minCnRegion({} -> {}) covered by hom-loss({})",
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

        int codingStart = dgData.TransData.CodingStart;
        int codingEnd = dgData.TransData.CodingEnd;

        for(final LohEvent lohEvent : mCopyNumberData.getLohData())
        {
            if(!lohEvent.Chromosome.equals(dgData.GeneData.Chromosome))
                continue;

            if(lohEvent.PosStart > codingEnd || lohEvent.PosEnd < codingStart)
                continue;

            LNX_LOGGER.debug("gene({}) coding region({} -> {}) covered by LOH({})",
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
                DriverEventType eventType = lohEvent.telomereLoss() ? LOH_SV_TELO : LOH_SV_CENTRO;
                DriverGeneEvent event = new DriverGeneEvent(eventType);
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

        LNX_LOGGER.debug("gene({}) chromosome({}) position({} -> {}) minCN({})",
                dgData.GeneData.GeneName, dgData.GeneData.Chromosome, dgData.TransData.TransStart, dgData.TransData.TransEnd,
                formatPloidy(dgData.GeneCN.minCopyNumber()));

        checkChromosomeAmplification(dgData);
        checkClusterAmplification(dgData, breakendList);

        if(dgData.getEvents().isEmpty())
        {
            LNX_LOGGER.debug("gene({}) AMP gain no event found", dgData.GeneData.GeneName);

            // otherwise no event
            DriverGeneEvent event = new DriverGeneEvent(GAIN);
            dgData.addEvent(event);
        }
    }

    private void checkChromosomeAmplification(final DriverGeneData dgData)
    {
        // check for an arm or whole chromosome amplified above the sample ploidy
        final TelomereCentromereCnData tcData = mCopyNumberData.getChrTeleCentroData().get(dgData.GeneData.Chromosome);

        if(tcData == null)
            return;

        double telomereCopyNumber = min(tcData.TelomerePArm, tcData.TelomereQArm);

        double centromereCopyNumber = min(tcData.CentromerePArm, tcData.CentromereQArm);

        double chromosomeCopyNumber = min(centromereCopyNumber, telomereCopyNumber);

        if(!isShortArmChromosome(dgData.GeneData.Chromosome))
        {
            if (chromosomeCopyNumber / mSamplePloidy > 2)
            {
                LNX_LOGGER.debug("gene({}) AMP gain from chromosome chChange({} telo={} centro={}) vs samplePloidy({})",
                        dgData.GeneData.GeneName, formatPloidy(chromosomeCopyNumber),
                        formatPloidy(telomereCopyNumber), formatPloidy(centromereCopyNumber), formatPloidy(mSamplePloidy));

                DriverGeneEvent event = new DriverGeneEvent(GAIN_CHR);
                event.setCopyNumberGain(chromosomeCopyNumber - mSamplePloidy);
                dgData.addEvent(event);
            }
        }

        double centromereCNChange = dgData.Arm == P_ARM ?
                min(tcData.CentromerePArm, tcData.TelomerePArm) - tcData.CentromereQArm
                : min(tcData.CentromereQArm, tcData.TelomereQArm) - tcData.CentromerePArm;

        // double centromereCNChange = dgData.Arm == CHROMOSOME_ARM_P ?
        //      tcData.CentromerePArm - tcData.CentromereQArm : tcData.CentromereQArm - tcData.CentromerePArm;

        if (centromereCNChange > 0 && !copyNumbersEqual(centromereCNChange, 0))
        {
            LNX_LOGGER.debug("gene({}) AMP gain from arm({}) cnChange across centromere({} -> {} = {})",
                    dgData.GeneData.GeneName, dgData.Arm, formatPloidy(tcData.CentromerePArm), formatPloidy(tcData.CentromereQArm),
                    formatPloidy(centromereCNChange));

            DriverGeneEvent event = new DriverGeneEvent(GAIN_ARM);
            event.setCopyNumberGain(centromereCNChange);
            dgData.addEvent(event);
        }
    }

    private OpposingSegment findOrCreateOpposingSegment(
            List<OpposingSegment> opposingSegments, final List<SvBreakend> breakendList, final SvBreakend startBreakend,
            boolean traverseUp, int stopPosition)
    {
        OpposingSegment opposingSegment = opposingSegments.stream()
                .filter(x -> x.Breakends.contains(startBreakend))
                .findFirst().orElse(null);

        if(opposingSegment != null)
            return opposingSegment;

        double startCopyNumber = startBreakend.getCopyNumber(traverseUp);
        double endCopyNumber = startBreakend.getCopyNumber(!traverseUp);;

        int index = startBreakend.getChrPosIndex();

        List<SvBreakend> segmentBreakends = Lists.newArrayList(startBreakend);

        while (true)
        {
            index += traverseUp ? 1 : -1;

            if (index < 0 || index >= breakendList.size())
                break;

            SvBreakend breakend = breakendList.get(index);

            if ((traverseUp && breakend.position() > stopPosition) || (!traverseUp && breakend.position() < stopPosition))
            {
                break;
            }

            if (breakend.getCluster() == startBreakend.getCluster())
            {
                segmentBreakends.add(breakend);
                endCopyNumber = breakend.getCopyNumber(!traverseUp);
            }
            else
            {
                break;
            }
        }

        double netClusterCNChange = endCopyNumber - startCopyNumber;

        if(netClusterCNChange > 0 || copyNumbersEqual(netClusterCNChange, 0))
        {
            opposingSegment = new OpposingSegment(startBreakend.getCluster(), segmentBreakends, 0);
            opposingSegments.add(opposingSegment);
            return null;
        }

        opposingSegment = new OpposingSegment(startBreakend.getCluster(), segmentBreakends, -netClusterCNChange);

        LNX_LOGGER.trace("added opposing segment: {}", opposingSegment);

        opposingSegments.add(opposingSegment);
        return opposingSegment;
    }

    private DriverAmpData checkClusterForAmplification(
            final DriverGeneData dgData, final List<SvBreakend> breakendList, final SvBreakend startBreakend, boolean traverseUp,
            List<OpposingSegment> opposingSegments)
    {
        int transStart = dgData.TransData.TransStart;
        int transEnd = dgData.TransData.TransEnd;

        double startCopyNumber = startBreakend.getCopyNumber(traverseUp);

        final SvCluster targetCluster = startBreakend.getCluster();

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

            if(index < 0 || index >= breakendList.size())
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
                    double endCopyNumber = 0;
                    if(breakend.arm().equals(segStartBreakend.arm()))
                    {
                        endCopyNumber = breakend.getCopyNumber(traverseUp);
                    }
                    else
                    {
                        // take the end copy number from the centromere
                        int prevIndex = index + (traverseUp ? -1 : 1);
                        final SvBreakend prevBreakend = breakendList.get(prevIndex);
                        endCopyNumber = prevBreakend.getCopyNumber(!traverseUp);
                    }

                    double segmentCNChange = endCopyNumber - segmentStartCopyNumber;
                    netClusterCNChange += segmentCNChange;

                    LNX_LOGGER.trace("gene({}) cluster({}) adding segment CN({} -> {} chg={}) net({}) breakends({} -> {})",
                            dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(segmentStartCopyNumber), formatPloidy(endCopyNumber),
                            formatPloidy(segmentCNChange), formatPloidy(netClusterCNChange), segStartBreakend, breakend);

                    inSegment = false;
                }

                OpposingSegment opposingSegment = findOrCreateOpposingSegment(
                        opposingSegments, breakendList, breakend, traverseUp, traverseUp ? transStart : transEnd);

                if(netClusterCNChange > 0 && opposingSegment != null && opposingSegment.remainingCNChange() > 0)
                {
                    if(opposingSegment.remainingCNChange() > netClusterCNChange)
                    {
                        LNX_LOGGER.trace("gene({}) cluster({}) netCN({}) cancelled by breakend({}) cnChange(orig={} remain={})",
                                dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(netClusterCNChange),
                                breakend, formatPloidy(breakend.copyNumberChange()), formatPloidy(opposingSegment.remainingCNChange()));

                        opposingSegment.reduceCNChange(netClusterCNChange);
                        netClusterCNChange = 0;
                    }
                    else
                    {
                        LNX_LOGGER.trace("gene({}) cluster({}) netCN({}) reducing by breakend({}) cnChange({})",
                                dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(netClusterCNChange),
                                breakend, formatPloidy(breakend.copyNumberChange()));

                        netClusterCNChange -= opposingSegment.remainingCNChange();
                        opposingSegment.zeroCNChange(); // keep so it won't be registered again but zero out
                    }
                }
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

            LNX_LOGGER.trace("gene({}) cluster({}) open segment startCN({}) net({}) start breakend({})",
                    dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(segmentStartCopyNumber), formatPloidy(clusterCNChange),
                    segStartBreakend);
        }

        if(clusterCNChange > 0 && !copyNumbersEqual(clusterCNChange, 0) && !copyNumbersEqual(startCopyNumber + clusterCNChange, startCopyNumber))
        {
            LNX_LOGGER.debug("gene({}) cluster({}) copy number gain({}) vs startCN({}) segments({}: CN={}) traversal({})",
                    dgData.GeneData.GeneName, targetCluster.id(), formatPloidy(clusterCNChange), formatPloidy(startCopyNumber),
                    segmentCount, formatPloidy(netClusterCNChange), traverseUp ? "up" : "down");

            return new DriverAmpData(targetCluster, traverseUp, breakendCount, segmentCount, startCopyNumber, clusterCNChange);
        }

        return null;
    }

    private static final double MIN_AMP_PERCENT_VS_MAX = 0.33;

    private void checkClusterAmplification(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        if (breakendList == null || breakendList.isEmpty())
            return;

        // sum up breakend ploidies from telomere to centromere for the gene in question net off ploidy within a cluster
        int transStart = dgData.TransData.TransStart;
        int transEnd = dgData.TransData.TransEnd;

        int startIndex = 0;
        int endIndex = breakendList.size() - 1;

        Map<SvCluster,DriverAmpData> clusterAmpData = Maps.newHashMap();

        for(int i = 0; i <= 1; ++i)
        {
            boolean traverseUp = (i == 0);

            List<SvCluster> processedClusters = Lists.newArrayList();
            List<OpposingSegment> opposingSegments = Lists.newArrayList();

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
                DriverAmpData ampData = checkClusterForAmplification(dgData, breakendList, breakend, traverseUp, opposingSegments);

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

        int index = 0;
        while(index < dgData.getEvents().size())
        {
            DriverGeneEvent event = dgData.getEvents().get(index);
            double armChrGain = event.getCopyNumberGain();

            if(armChrGain < MIN_AMP_PERCENT_VS_MAX * maxCnChange)
            {
                dgData.getEvents().remove(index);
            }
            else
            {
                maxCnChange = max(maxCnChange, armChrGain);
                ++index;
            }
        }

        for(Map.Entry<SvCluster,DriverAmpData> entry : clusterAmpData.entrySet())
        {
            final SvCluster cluster = entry.getKey();
            final DriverAmpData ampData = entry.getValue();

            // take any at least 20% of the largest contributing cluster
            if(ampData.NetCNChange / maxCnChange < MIN_AMP_PERCENT_VS_MAX)
                continue;

            LNX_LOGGER.debug("gene({}) cluster({}) adding AMP data: {}",
                    dgData.GeneData.GeneName, cluster.id(), ampData);

            DriverGeneEvent event = new DriverGeneEvent(GAIN);
            event.setCluster(cluster);
            event.setCopyNumberGain(ampData.NetCNChange);
            event.setAmpData(ampData);

            event.setSvInfo(String.format("%s;%d;%d;%.1f",
                    ampData.TraverseUp, ampData.BreakendCount, ampData.SegmentCount, ampData.StartCopyNumber));
            dgData.addEvent(event);
        }
    }

    private void findDisruptiveDelDrivers()
    {
        final List<String> delDriverGeneIds = mDriverGeneDataList.stream()
                .filter(x -> x.DriverData.driver() == DriverType.DEL)
                .map(x -> x.GeneData.GeneId).collect(Collectors.toList());

        List<DriverGeneData> disDelDrivers = Lists.newArrayList();

        for(Map.Entry<String,List<SvBreakend>> entry : mChrBreakendMap.entrySet())
        {
            for(final SvBreakend breakend : entry.getValue())
            {
                boolean delType = breakend.orientation() == 1 && breakend.getDBLink() != null;
                boolean dupType = breakend.orientation() == -1 && breakend.getSV().type() == DUP;

                if(!delType && !dupType)
                    continue;

                List<GeneAnnotation> genesList = breakend.getSV().getGenesList(breakend.usesStart()).stream()
                        .filter(x -> !delDriverGeneIds.contains(x.StableId))
                        .filter(x -> mReportableDelGeneIds.contains(x.StableId))
                        .collect(Collectors.toList());

                if(genesList.isEmpty())
                    continue;

                // consider any disruptive canonical transcript
                for(final GeneAnnotation gene : genesList)
                {
                    final Transcript trans = gene.canonical();

                    if(trans == null || !trans.isDisruptive())
                        continue;

                    if(delType)
                    {
                        final SvLinkedPair dbLink = breakend.getDBLink();

                        // calculate the copy number over the deletion bridge section
                        double cnLowSide = breakend.copyNumberLowSide();

                        double otherSvPloidy = dbLink.getOtherBreakend(breakend).ploidy();

                        // account for an overlapping DB by subtracting the ploidy of the overlap
                        if(dbLink.length() < 0)
                            cnLowSide -= otherSvPloidy;

                        if(cnLowSide < TOTAL_CN_LOSS)
                        {
                            delDriverGeneIds.add(gene.StableId);

                            LNX_LOGGER.debug("gene({}) cluster({}) breakend({}) cause homozyous disruption for cnLowSide({}) dbLength({}) otherSvPloidy({})",
                                    trans.geneName(), breakend.getCluster().id(), breakend,
                                    formatPloidy(cnLowSide), dbLink.length(), formatPloidy(otherSvPloidy));

                            DriverGeneData dgData = createDriverData(gene);
                            DriverGeneEvent event = new DriverGeneEvent(HOM_DEL_DISRUPTION);
                            event.addSvBreakendPair(dbLink.firstBreakend(), dbLink.secondBreakend(), "DB");
                            event.setCluster(breakend.getCluster());
                            dgData.addEvent(event);
                            disDelDrivers.add(dgData);
                        }
                    }
                    else
                    {
                        // DUP must be wholy contained within the same gene and
                        if(!breakend.getSV().getGenesList(!breakend.usesStart()).stream().anyMatch(x -> x.StableId == gene.StableId))
                            continue;

                        final SvBreakend otherBreakend = breakend.getOtherBreakend();

                        double cnLowSideStart = breakend.copyNumberLowSide();
                        double cnLowSideEnd = otherBreakend.copyNumberLowSide();
                        double ploidy = breakend.ploidy();
                        double ploidyThreshold = max(ploidy * (1 + MAX_COPY_NUM_DIFF_PERC), ploidy + MAX_COPY_NUM_DIFF);

                        if(cnLowSideStart < ploidyThreshold && cnLowSideEnd < ploidyThreshold)
                        {
                            delDriverGeneIds.add(gene.StableId);

                            LNX_LOGGER.debug("gene({}) cluster({}) DUP({}) cause homozygous disruption cnLowSide({} & {}) ploidy({})",
                                    trans.geneName(), breakend.getCluster().id(), breakend.getSV().id(),
                                    formatPloidy(cnLowSideStart), formatPloidy(cnLowSideEnd), formatPloidy(ploidy));

                            DriverGeneData dgData = createDriverData(gene);
                            DriverGeneEvent event = new DriverGeneEvent(HOM_DUP_DISRUPTION);
                            event.addSvBreakendPair(breakend, otherBreakend, "DUP");
                            event.setCluster(breakend.getCluster());
                            dgData.addEvent(event);
                            disDelDrivers.add(dgData);
                        }
                    }
                }
            }
        }

        disDelDrivers.forEach(x -> writeDriverData(x));
    }

    private DriverGeneData createDriverData(final GeneAnnotation gene)
    {
        GeneCopyNumber gcnData = null;

        if(mDbAccess != null)
        {
            final List<GeneCopyNumber> geneCopyNumbers = mDbAccess.readGeneCopynumbers(mSampleId, Lists.newArrayList(gene.GeneName));
            gcnData = !geneCopyNumbers.isEmpty() ? geneCopyNumbers.get(0) : null;
        }
        else
        {
            gcnData = findGeneCopyNumber(gene.GeneName);
        }

        final DriverCatalog driverRecord = ImmutableDriverCatalog.builder()
                .driver(DriverType.HOM_DISRUPTION)
                .category(TSG)
                .gene(gene.GeneName)
                .chromosome(gene.chromosome())
                .chromosomeBand(gene.karyotypeBand())
                .likelihoodMethod(DEL)
                .driverLikelihood(1.0)
                .dndsLikelihood(0)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(false)
                .minCopyNumber(gcnData != null ? gcnData.minCopyNumber() : 0)
                .maxCopyNumber(gcnData != null ? gcnData.maxCopyNumber() : 0)
                .build();

        mDriverCatalog.add(driverRecord);

        final EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(gene.GeneName);

        if (geneData == null)
            return null;

        final TranscriptData canonicalTrans = mGeneTransCache.getTranscriptData(geneData.GeneId, "");

        DriverGeneData dgData = new DriverGeneData(driverRecord, geneData, canonicalTrans, gcnData);
        mDriverGeneDataList.add(dgData);

        return dgData;
    }

    private final GeneCopyNumber findGeneCopyNumber(final String geneName)
    {
        for(GeneCopyNumber geneCN : mGeneCopyNumberData)
        {
            if(geneCN.gene().equals(geneName))
                return geneCN;
        }

        return null;
    }

    private void writeSampleDriverData()
    {
        if(mConfig.hasMultipleSamples())
            return;

        if(mDriverOutputList.isEmpty())
            return;

        if(mConfig.UploadToDB && mDbAccess != null)
        {
            mDbAccess.writeDriverCatalog(mSampleId, mDriverCatalog);

            mDbAccess.writeSvDrivers(mSampleId, mDriverOutputList);
        }

        try
        {
            final String driverCatalogFile = DriverCatalogFile.generateFilename(mOutputDir, mSampleId);
            DriverCatalogFile.write(driverCatalogFile, mDriverCatalog);

            final String driversFile = LinxDriverFile.generateFilename(mOutputDir, mSampleId);
            LinxDriverFile.write(driversFile, mDriverOutputList);
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to write drivers file: {}", e.toString());
        }
    }

    private void writeDriverData(final DriverGeneData dgData)
    {
        // convert to a sample driver record
        if(mVisWriter != null)
        {
            for (final DriverGeneEvent driverEvent : dgData.getEvents())
            {
                int clusterId = driverEvent.getCluster() != null ? driverEvent.getCluster().id() : -1;

                mDriverOutputList.add(ImmutableLinxDriver.builder()
                        .clusterId(clusterId)
                        .gene(dgData.GeneData.GeneName)
                        .eventType(driverEvent.Type.toString())
                        .build());

                mVisWriter.addGeneExonData(clusterId, dgData.GeneData.GeneId, dgData.GeneData.GeneName,
                        "", 0, dgData.GeneData.Chromosome, GENE_TYPE_DRIVER);
            }
        }

        if(!mConfig.hasMultipleSamples() || mOutputDir.isEmpty())
            return;

        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += "LNX_DRIVERS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Gene,Category,DriverType,LikelihoodMethod,Likelihood");
                mFileWriter.write(",FullyMatched,EventType,ClusterId,ClusterCount,ResolvedType");
                mFileWriter.write(",Chromosome,Arm,SamplePloidy,GeneMinCN,CentromereCN,TelomereCN,CNGain");
                mFileWriter.write(",SvIdStart,SvIdEnd,SvPosStart,SvPosEnd,SvMatchType");
                mFileWriter.newLine();
            }

            BufferedWriter writer = mFileWriter;

            final DriverCatalog driverGene = dgData.DriverData;
            final EnsemblGeneData geneData = dgData.GeneData;

            final TelomereCentromereCnData tcData = mCopyNumberData.getChrTeleCentroData().get(dgData.GeneData.Chromosome);

            double centromereCopyNumber = dgData.Arm == P_ARM ? tcData.CentromerePArm : tcData.CentromereQArm;
            double telomereCopyNumber = dgData.Arm == P_ARM ? tcData.TelomerePArm : tcData.TelomereQArm;

            for(final DriverGeneEvent driverEvent : dgData.getEvents())
            {
                final SvBreakend[] breakendPair = driverEvent.getBreakendPair();

                writer.write(String.format("%s,%s,%s,%s,%s,%.2f,%s,%s",
                        mSampleId, driverGene.gene(), driverGene.category(), driverGene.driver(),
                        driverGene.likelihoodMethod(), driverGene.driverLikelihood(), dgData.fullyMatched(), driverEvent.Type));

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

                // Chromosome,Arm,SamplePloidy,GeneMinCN,CentromereCN,TelomereCN,CNGain
                // SvIdStart,SvIdEnd,SvPosStart,SvPosEnd,SvMatchType

                writer.write(String.format(",%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f",
                        geneData.Chromosome, dgData.Arm, mSamplePloidy, dgData.GeneCN != null ? dgData.GeneCN.minCopyNumber() : -1,
                        centromereCopyNumber, telomereCopyNumber, driverEvent.getCopyNumberGain()));

                int posStart = breakendPair[SE_START] != null ? breakendPair[SE_START].position() : 0;
                int posEnd = breakendPair[SE_END] != null ? breakendPair[SE_END].position() : 0;
                String svIdStart = breakendPair[SE_START] != null ? breakendPair[SE_START].getSV().idStr() : "-1";
                String svIdEnd = breakendPair[SE_END] != null ? breakendPair[SE_END].getSV().idStr() : "-1";

                writer.write(String.format(",%s,%s,%d,%d,%s",
                        svIdStart, svIdEnd, posStart, posEnd, driverEvent.getSvInfo()));


                writer.newLine();
            }
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing driver data to outputFile: {}", e.toString());
        }
    }

    @VisibleForTesting
    public void addDriverGene(final DriverCatalog driver, final GeneCopyNumber geneCopyNumber)
    {
        mDriverCatalog.clear();
        mDriverCatalog.add(driver);

        mGeneCopyNumberData.clear();
        mGeneCopyNumberData.add(geneCopyNumber);
    }

    public void close()
    {
        if(LNX_LOGGER.isDebugEnabled() || mConfig.hasMultipleSamples())
        {
            mPerfCounter.logStats();
        }

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
                LNX_LOGGER.error("Empty gene copy number CSV file({})", gcnFileName);
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
                    LNX_LOGGER.error("GCN file invalid item count({}) vs expected({})", items.length, GENE_CN_DATA_FILE_ITEM_COUNT);
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

            LNX_LOGGER.info("loaded {} gene copy-number records from file", rowCount);
        }
        catch (IOException e)
        {
            LNX_LOGGER.error("Failed to read gene copy number CSV file({}): {}", gcnFileName, e.toString());
        }
    }

}

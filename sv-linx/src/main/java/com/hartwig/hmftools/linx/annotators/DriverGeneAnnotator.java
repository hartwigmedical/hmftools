package com.hartwig.hmftools.linx.annotators;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArm;
import static com.hartwig.hmftools.linx.annotators.DriverEventType.GAIN;
import static com.hartwig.hmftools.linx.annotators.DriverEventType.GAIN_ARM;
import static com.hartwig.hmftools.linx.annotators.DriverEventType.GAIN_CHR;
import static com.hartwig.hmftools.linx.annotators.DriverEventType.LOH;
import static com.hartwig.hmftools.linx.annotators.DriverEventType.LOH_ARM;
import static com.hartwig.hmftools.linx.annotators.DriverEventType.LOH_CHR;
import static com.hartwig.hmftools.linx.annotators.DriverGeneEvent.SV_DRIVER_TYPE_ARM_SV;
import static com.hartwig.hmftools.linx.annotators.DriverGeneEvent.SV_DRIVER_TYPE_DEL;
import static com.hartwig.hmftools.linx.annotators.DriverGeneEvent.SV_DRIVER_TYPE_DUP;
import static com.hartwig.hmftools.linx.annotators.DriverGeneEvent.SV_DRIVER_TYPE_FB;
import static com.hartwig.hmftools.linx.annotators.DriverGeneEvent.SV_DRIVER_TYPE_TI;
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
import java.sql.Timestamp;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.ImmutableGeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.cn.LohEvent;
import com.hartwig.hmftools.linx.types.SvLinkedPair;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class DriverGeneAnnotator
{
    private static final Logger LOGGER = LogManager.getLogger(DriverGeneAnnotator.class);

    final DatabaseAccess mDbAccess;
    SvGeneTranscriptCollection mGeneTransCache;

    private List<DriverCatalog> mDriverCatalog;
    private List<GeneCopyNumber> mGeneCopyNumberData;
    private BufferedWriter mFileWriter;
    private String mOutputDir;
    private double mSamplePloidy;
    private List<DriverGeneData> mDriverGeneDataList;

    private PerformanceCounter mPerfCounter;

    // references only
    private String mSampleId;
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private List<LohEvent> mLohEventList;
    private List<HomLossEvent> mHomLossList;
    private Map<String, double[]> mChrCopyNumberMap;
    private Map<String, List<GeneCopyNumber>> mSampleGeneCopyNumberMap;
    private VisualiserWriter mVisWriter;

    // temp
    private boolean mWriteMatchedGeneCopyNumber;
    private BufferedWriter mGCNFileWriter;

    public DriverGeneAnnotator(DatabaseAccess dbAccess, SvGeneTranscriptCollection geneTranscriptCollection, final String outputDir)
    {
        mDbAccess = dbAccess;
        mGeneTransCache = geneTranscriptCollection;

        mDriverCatalog = Lists.newArrayList();
        mGeneCopyNumberData = Lists.newArrayList();
        mDriverGeneDataList = Lists.newArrayList();
        mSampleGeneCopyNumberMap = new HashMap();
        mChrCopyNumberMap = null;
        mFileWriter = null;
        mOutputDir = outputDir;
        mSamplePloidy = 0;

        mGeneTransCache.createGeneNameIdMap();

        mWriteMatchedGeneCopyNumber = false;
        mGCNFileWriter = null;
        mVisWriter = null;

        mPerfCounter = new PerformanceCounter("Drivers");
    }

    private static final String WRITE_GCN_DATA = "write_gcn_data";
    private static final String GCN_DATA_FILE = "gcn_data_file";

    public static void addCmdLineArgs(Options options)
    {
        options.addOption(WRITE_GCN_DATA, false, "Write a cache of driver-matched gene copy number data");
        options.addOption(GCN_DATA_FILE, true, "Cache of driver-matched gene copy number data");
    }

    public boolean loadConfig(final CommandLine cmd)
    {
        mWriteMatchedGeneCopyNumber = cmd.hasOption(WRITE_GCN_DATA);

        initialiseGeneData(cmd.getOptionValue(GCN_DATA_FILE,""));

        return true;
    }

    public void setCopyNumberData(
            final Map<String, double[]> chrCopyNumberMap,
            final List<LohEvent> lohData,
            final List<HomLossEvent> homLossData)
    {
        mLohEventList = lohData;
        mHomLossList = homLossData;
        mChrCopyNumberMap = chrCopyNumberMap;
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

        LOGGER.debug("sample({}) retrieved {} driver gene records", sampleId, mDriverCatalog.size());
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

        if(!mWriteMatchedGeneCopyNumber)
        {
            mGeneCopyNumberData = mDbAccess.readGeneCopynumbers(sampleId);
            return;
        }

        LOGGER.info("sample({}) writing gene copy number data", mSampleId);

        List<GeneCopyNumber> gcnDataList = mDbAccess.readGeneCopynumbers(sampleId);

        for(final DriverCatalog driverData : mDriverCatalog)
        {
            for(final GeneCopyNumber gcnData : gcnDataList)
            {
                if(gcnData.gene().equals(driverData.gene()))
                {
                    writeGeneCopyNumberData(gcnData);
                    mGeneCopyNumberData.add(gcnData);
                    break;
                }
            }
        }

        LOGGER.info("sample({}) drivers matched({}/{}) from {} gene copy number data",
                mSampleId, mGeneCopyNumberData.size(), mDriverCatalog.size(), gcnDataList.size());
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

            final GeneCopyNumber geneCN = findGeneCopyNumber(driverGene);

            if (geneCN == null)
            {
                LOGGER.warn("sample({}) gene({}) copy number data not found", mSampleId, driverGene.gene());
                continue;
            }

            DriverGeneData dgData = new DriverGeneData(driverGene, geneData, canonicalTrans, geneCN);
            mDriverGeneDataList.add(dgData);

            if(dgData.DriverData.category() == TSG && dgData.DriverData.driver() == DriverType.DEL)
            {
                if(dgData.DriverData.likelihoodMethod() == LikelihoodMethod.DEL)
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

                if (breakendList == null || breakendList.isEmpty())
                {
                    writeDriverData(dgData);
                    continue;
                }

                annotateAmplification(dgData, breakendList);
                // annotateSpanningSVs(dgData, breakendList);
            }
        }

        mChrBreakendMap = null;

        mPerfCounter.stop();
    }

    private boolean driverTypeHandled(final DriverCatalog driverGene)
    {
        if(driverGene.category() == TSG && driverGene.driver() == DriverType.DEL)
            return true;

        if(driverGene.category() == TSG && driverGene.likelihoodMethod() == LikelihoodMethod.BIALLELIC)
            return true;

        if(driverGene.category() == ONCO && driverGene.driver() == DriverType.AMP)
            return true;

        return false;
    }

    private void annotateDeleteEvent(final DriverGeneData dgData)
    {
        // find the LOH and hom-loss events which caused this DEL
        long minRegionStart = dgData.GeneCN.minRegionStart();
        long minRegionEnd = dgData.GeneCN.minRegionEnd();

        for(final LohEvent lohEvent : mLohEventList)
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

        writeDriverData(dgData);
    }

    private void annotateBiallelicEvent(DriverGeneData dgData)
    {
        // look for an LOH covering the min-gene region
        long minRegionStart = dgData.GeneCN.minRegionStart();
        long minRegionEnd = dgData.GeneCN.minRegionEnd();

        for(final LohEvent lohEvent : mLohEventList)
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
                    DriverGeneEvent event = new DriverGeneEvent(LOH_ARM);
                    event.setLohEvent(lohEvent);
                    event.addSvBreakendPair(lohEvent.getBreakend(true), lohEvent.getBreakend(false), SV_DRIVER_TYPE_ARM_SV);
                    dgData.addEvent(event);
                }
            }

            break;
        }

        writeDriverData(dgData);
    }

    private void annotateAmplification(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        annotateArmAmplification(dgData);
        annotateBreakendAmplification(dgData, breakendList);
    }

    private void annotateArmAmplification(final DriverGeneData dgData)
    {
        // check for an arm or whole chromosome amplified above the sample ploidy
        if(mSamplePloidy == 0)
            return;

        double[] cnData = mChrCopyNumberMap.get(dgData.GeneData.Chromosome);

        double telomereCopyNumber = min(cnData[P_ARM_TELOMERE_CN], cnData[Q_ARM_TELOMERE_CN]);

        double armCopyNumber = min(cnData[CENTROMERE_CN],
                dgData.Arm == CHROMOSOME_ARM_P ? cnData[P_ARM_TELOMERE_CN] : cnData[Q_ARM_TELOMERE_CN]);

        double chromosomeCopyNumber = min(cnData[CENTROMERE_CN], telomereCopyNumber);

        if(chromosomeCopyNumber/mSamplePloidy > 2)
        {
            dgData.addEvent(new DriverGeneEvent(GAIN_CHR));
        }
        else if(armCopyNumber/mSamplePloidy > 2)
        {
            dgData.addEvent(new DriverGeneEvent(GAIN_ARM));
        }
    }

    private void annotateBreakendAmplification(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        // find the cause - DUP, foldback, otherwise assume whole-chromatid duplication

        // trying to find the breakends which amplify this gene
        // take any INV, DUP or TI which straddles the gene

        final DriverCatalog driverGene = dgData.DriverData;
        long transStart = dgData.TransData.TransStart;
        long transEnd = dgData.TransData.TransEnd;

        double maxCopyNumber = 0;
        SvVarData maxSvStart = null;
        SvVarData maxSvEnd = null;
        SvBreakend foldbackBreakend = null;
        List<Integer> reportedClusters = Lists.newArrayList();

        for (int i = 0; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);

            if(breakend.position() > transStart)
            {
                // past the gene but AMP may be explained by an unchained foldback
                if(foldbackBreakend == null && breakend.orientation() == 1
                && !breakend.getSV().getFoldbackLink(breakend.usesStart()).isEmpty())
                {
                    foldbackBreakend = breakend;
                    break;
                }

                continue; // looking for a foldback as an explanation
            }

            if(breakend.orientation() == 1)
                continue;

            final SvVarData varStart = breakend.getSV();

            if(breakend.orientation() == -1 && !varStart.getFoldbackLink(breakend.usesStart()).isEmpty())
                foldbackBreakend = breakend;

            if(varStart.type() == DUP)
            {
                if(varStart.position(true) > transStart || varStart.position(false) < transEnd)
                    continue;

                // require the DUP to no be interupted by any other SV within this cluster within this gene region
                SvBreakend lowerBreakend = breakend;
                SvBreakend upperBreakend = varStart.getBreakend(false);

                if(lowerBreakend.getClusterChrPosIndex() != upperBreakend.getClusterChrPosIndex() - 1)
                    continue;

                LOGGER.debug(String.format("sample(%s) cluster(%s) gene(%s) single SV(%s %s) cn(%.2f) cnChg(%.2f)",
                        mSampleId, varStart.getCluster().id(), driverGene, varStart.posId(), varStart.type(),
                        varStart.copyNumber(true), varStart.copyNumberChange(true)));

                DriverGeneEvent event = new DriverGeneEvent(GAIN);
                event.addSvBreakendPair(lowerBreakend, upperBreakend, SV_DRIVER_TYPE_DUP);
                dgData.addEvent(event);

                if(varStart.copyNumberChange(true) > maxCopyNumber)
                {
                    maxCopyNumber = varStart.copyNumberChange(true);
                    maxSvStart = varStart;
                    maxSvEnd = varStart;
                }

                continue;
            }

            // only report on each cluster once
            if(reportedClusters.contains(varStart.getCluster().id()))
                continue;

            // look for the first TI which overlaps the gene region
            List<SvLinkedPair> tiPairs = varStart.getLinkedPairs(breakend.usesStart());

            for(SvLinkedPair tiPair : tiPairs)
            {
                if (breakend.position() + tiPair.length() < transEnd)
                    continue;

                final SvVarData varEnd = tiPair.getOtherSV(varStart);
                final SvBreakend beStart = tiPair.getBreakend(true).getOrigBreakend();
                final SvBreakend beEnd = tiPair.getBreakend(false).getOrigBreakend();

                LOGGER.debug(String.format("sample(%s) cluster(%d fb=%s) gene(%s) SVs start(%s cn=%.2f cnChg=%.2f) end(%s cn=%.2f cnChg=%.2f) in linked pair",
                        mSampleId, varStart.getCluster().id(), varStart.getCluster().getFoldbacks().size(), driverGene,
                        varStart.posId(), varStart.copyNumber(beStart.usesStart()), varStart.copyNumberChange(beStart.usesStart()),
                        varEnd.posId(), varEnd.copyNumber(beEnd.usesStart()), varEnd.copyNumberChange(beEnd.usesStart())));

                DriverGeneEvent event = new DriverGeneEvent(GAIN);
                event.addSvBreakendPair(beStart, beEnd, SV_DRIVER_TYPE_TI);
                dgData.addEvent(event);

                reportedClusters.add(varStart.getCluster().id());

                if (varStart.copyNumberChange(true) > maxCopyNumber)
                {
                    maxCopyNumber = varStart.copyNumberChange(true);
                    maxSvStart = varStart;
                    maxSvEnd = varEnd;
                }

                break;
            }
        }

        if(maxSvStart == null || maxSvEnd == null)
        {
            if(foldbackBreakend != null)
            {
                SvVarData foldbackSv = foldbackBreakend.getSV();
                SvBreakend otherBreakend = foldbackSv.isChainedFoldback() ?
                        foldbackSv.getChainedFoldbackSv().getChainedFoldbackBreakend() : foldbackBreakend.getOtherBreakend();

                if(otherBreakend != null)
                {
                    DriverGeneEvent event = new DriverGeneEvent(GAIN);
                    event.addSvBreakendPair(foldbackBreakend, otherBreakend, SV_DRIVER_TYPE_FB);
                    dgData.addEvent(event);
                }
            }
        }

        writeDriverData(dgData);
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

    private void writeDriverData(final DriverGeneData dgData)
    {
        try
        {
            if(mFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += "SVA_DRIVERS.csv";

                mFileWriter = createBufferedWriter(outputFileName, false);

                mFileWriter.write("SampleId,Gene,Category,DriverType,LikelihoodMethod");
                mFileWriter.write(",FullyMatched,EventType,ClusterId,ClusterCount,ResolvedType");
                mFileWriter.write(",PosStart,PosEnd,SvIdStart,SvIdEnd,MatchInfo");
                mFileWriter.write(",SamplePloidy,Chromosome,Arm,MinCN,CentromereCN,TelomereCN");
                mFileWriter.newLine();
            }

            BufferedWriter writer = mFileWriter;

            final DriverCatalog driverGene = dgData.DriverData;
            final EnsemblGeneData geneData = dgData.GeneData;
            double[] cnData = mChrCopyNumberMap.get(geneData.Chromosome);
            int refClusterId = -1;

            for(final DriverGeneEvent driverEvent : dgData.getEvents())
            {
                final SvBreakend[] breakendPair = driverEvent.getBreakendPair();

                writer.write(String.format("%s,%s,%s,%s,%s,%s,%s",
                        mSampleId, driverGene.gene(), driverGene.category(), driverGene.driver(), driverGene.likelihoodMethod(),
                        dgData.fullyMatched(), driverEvent.Type));

                final SvCluster cluster = driverEvent.getCluster();

                if(cluster != null)
                {
                    refClusterId = cluster.id();
                    writer.write(String.format(",%d,%d,%s", cluster.id(), cluster.getSvCount(), cluster.getResolvedType()));
                }
                else
                {
                    writer.write(String.format(",-1,0,"));
                }

                // breakend info if present

                long posStart = breakendPair[SE_START] != null ? breakendPair[SE_START].position() : 0;
                long posEnd = breakendPair[SE_END] != null ? breakendPair[SE_END].position() : 0;
                String svIdStart = breakendPair[SE_START] != null ? breakendPair[SE_START].getSV().id() : "-1";
                String svIdEnd = breakendPair[SE_END] != null ? breakendPair[SE_END].getSV().id() : "-1";

                writer.write(String.format(",%d,%d,%s,%s,%s",
                        posStart, posEnd, svIdStart, svIdEnd, driverEvent.getSvInfo()));


                writer.write(String.format(",%.2f,%s,%s,%.2f,%.2f,%.2f",
                        mSamplePloidy, geneData.Chromosome, dgData.Arm, dgData.GeneCN != null ? dgData.GeneCN.minCopyNumber() : -1,
                        cnData[CENTROMERE_CN], dgData.Arm == CHROMOSOME_ARM_P ? cnData[P_ARM_TELOMERE_CN] : cnData[Q_ARM_TELOMERE_CN]));

                writer.newLine();
            }

            mVisWriter.addGeneExonData(refClusterId, geneData.GeneId, geneData.GeneName,
                    "", 0, geneData.Chromosome, "DRIVER");

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
        closeBufferedWriter(mGCNFileWriter);
    }

    private static String GENE_CN_DATA_FILE = "SVA_GENE_COPY_NUMBER.csv";
    private static int GENE_CN_DATA_FILE_ITEM_COUNT = 30;

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

            while ((line = fileReader.readLine()) != null)
            {
                String[] items = line.split(",");

                if (items.length != GENE_CN_DATA_FILE_ITEM_COUNT)
                    continue;

                String sampleId = items[1];

                if(currentSample.isEmpty() || !currentSample.equals(sampleId))
                {
                    gcnDataList = Lists.newArrayList();
                    currentSample = sampleId;
                    mSampleGeneCopyNumberMap.put(currentSample, gcnDataList);

                }

                //                          0       1       2           3   4   5      6              7           8
//                mGCNFileWriter.write("modified,sampleId,chromosome,start,end,gene,chromosomeBand,transcriptId,transcriptVersion");
                //                       10             11           12             13                  14                 15        16              17
//                mGCNFileWriter.write(",minCopyNumber,maxCopyNumber,somaticRegions,germlineHomRegions,germlineHetRegions,minRegions,minRegionStart,minRegionEnd");
                //                          18                  19                      20          21
//                mGCNFileWriter.write(",minRegionStartSupport,minRegionEndSupport,minRegionMethod,minMinorAllelePloidy");

                int index = 2;

                GeneCopyNumber gcnData = ImmutableGeneCopyNumber.builder()
                        .chromosome(items[index++])
                        .start(Integer.parseInt(items[index++]))
                        .end(Integer.parseInt(items[index++]))
                        .gene(items[index++])
                        .chromosomeBand(items[index++])
                        .transcriptID(items[index++])
                        .transcriptVersion(Integer.parseInt(items[index++]))
                        .minCopyNumber(Double.parseDouble(items[index++]))
                        .maxCopyNumber(Double.parseDouble(items[index++]))
                        .somaticRegions(Integer.parseInt(items[index++]))
                        .germlineHomRegions(Integer.parseInt(items[index++]))
                        .germlineHet2HomRegions(Integer.parseInt(items[index++]))
                        .minRegions(Integer.parseInt(items[index++]))
                        .minRegionStart(Integer.parseInt(items[index++]))
                        .minRegionEnd(Integer.parseInt(items[index++]))
                        .minRegionStartSupport(SegmentSupport.valueOf(items[index++]))
                        .minRegionEndSupport(SegmentSupport.valueOf(items[index++]))
                        .minRegionMethod(CopyNumberMethod.valueOf(items[index++]))
                        .minMinorAllelePloidy(Double.parseDouble(items[index++]))
                        .build();


                gcnDataList.add(gcnData);
            }
        }
        catch (IOException e)
        {
            LOGGER.error("Failed to read gene copy number CSV file({}): {}", gcnFileName, e.toString());
        }
    }

    private void writeGeneCopyNumberData(final GeneCopyNumber gcnData)
    {
        try
        {
            if(mGCNFileWriter == null)
            {
                String outputFileName = mOutputDir;

                outputFileName += GENE_CN_DATA_FILE;

                mGCNFileWriter = createBufferedWriter(outputFileName, false);

                mGCNFileWriter.write("modified,sampleId,chromosome,start,end,gene,chromosomeBand,transcriptId,transcriptVersion");
                mGCNFileWriter.write(",minCopyNumber,maxCopyNumber,somaticRegions,germlineHomRegions,germlineHetRegions,minRegions,minRegionStart,minRegionEnd");
                mGCNFileWriter.write(",minRegionStartSupport,minRegionEndSupport,minRegionMethod,minMinorAllelePloidy");
                mGCNFileWriter.newLine();
            }

            BufferedWriter writer = mGCNFileWriter;

            final Timestamp timestamp = new Timestamp(new Date().getTime());

            writer.write(String.format("%s,%s,%s,%d,%d,%s,%s,%s,%d",
                    timestamp, mSampleId, gcnData.chromosome(), gcnData.start(), gcnData.end(),
                    gcnData.gene(), gcnData.chromosomeBand(), gcnData.transcriptID(), gcnData.transcriptVersion()));

            writer.write(String.format(",%.4f,%.4f,%d,%d,%d,%d,%d,%d,%s,%s,%s",
                    gcnData.minCopyNumber(), gcnData.maxCopyNumber(), gcnData.somaticRegions(),
                    gcnData.germlineHomRegions(), gcnData.germlineHet2HomRegions(),
                    gcnData.minRegions(), gcnData.minRegionStart(), gcnData.minRegionEnd(),
                    gcnData.minRegionStartSupport(), gcnData.minRegionEndSupport(), gcnData.minRegionMethod().toString()));


            writer.write(String.format(",%.4f,", gcnData.minMinorAllelePloidy()));

            writer.newLine();

        }
        catch (final IOException e)
        {
            LOGGER.error("error writing gene copy number to outputFile: {}", e.toString());
        }
    }

    private void annotateDeleteEvent_old(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        /* DEL identification:
            - 1 or 2 SVs which caused this, start from DEL region (ie gene) and work out
            - gene copy number table - minCopyRegion < 0.5 makes it a DEL-type driver candidate
            - walk out from here to end of LOH event ie where (1- BAF*CN) > 0.5 again
            - see COLO829T and PTEN as an example - here there is a simple DEL and then most of the other chromatid has been lost
         */

    /*
        // first find the min copy number within the gene region
        // then walk out in both directions to find the SV which caused the loss
        // and then walk out again until heterozygosity is gained (ie the end of an LOH)

        SvBreakend minBreakend = null; // breakend with lowest copy number covering any part of the gene region

        final DriverCatalog driverGene = dgData.DriverData;
        final GeneCopyNumber geneCN = dgData.GeneCN;

        for (int i = 0; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);
            final SvBreakend nextBreakend = i < breakendList.size() - 1 ? breakendList.get(i + 1) : null;

            if(breakend.orientation() == -1)
            {
                if(breakend.position() >= geneCN.minRegionEnd())
                {
                    if(minBreakend == null)
                    {
                        // first breakend is already past the gene's min copy number region
                        minBreakend = breakend;
                    }

                    break;
                }

                continue;
            }

            if(breakend.position() >= geneCN.minRegionEnd())
            {
                LOGGER.debug("sample({}) gene({}) breakend past min-gene region({} - {})",
                        mSampleId, driverGene, geneCN.minRegionStart(), geneCN.minRegionEnd());
                break;
            }

            // not at the correct region yet
            if(breakend.position() < geneCN.minRegionStart() && nextBreakend != null && nextBreakend.position() < geneCN.minRegionStart())
                continue;

            // find the breakend closest to the min region with positive orientation
            if(breakend.position() <= geneCN.minRegionStart() && nextBreakend != null && nextBreakend.position() > geneCN.minRegionStart())
            {
                minBreakend = breakend;
            }
            else if(breakend.position() <= geneCN.minRegionStart() && nextBreakend == null)
            {
                minBreakend = breakend;
            }
            else if(minBreakend != null && breakend.position() < geneCN.minRegionEnd())
            {
                // a tighter breakend to the one allocated
                minBreakend = breakend;
            }
        }

        if(minBreakend == null)
        {
            LOGGER.debug("sample({}) gene({}) not allocated to SVs", mSampleId, driverGene);
            writeDriverData(dgData);
            return;
        }

        LOGGER.debug(String.format("gene(%s) min copy number at breakend(%s cn=%.2f cnChg=%.2f)",
                driverGene.gene(), minBreakend.toString(), minBreakend.getCopyNumber(false),
                minBreakend.getSV().copyNumberChange(minBreakend.usesStart())));

        SvBreakend startBreakend = minBreakend.orientation() == 1 ? minBreakend : null;
        SvBreakend endBreakend = null;

        if(startBreakend != null)
        {
            for(int index = startBreakend.getChrPosIndex() + 1; index < breakendList.size(); ++index)
            {
                final SvBreakend breakend = breakendList.get(index);

                if(breakend.orientation() == -1 && breakend.getCluster() == startBreakend.getCluster())
                {
                    endBreakend = breakend;
                    break;
                }
            }

            // now walk forwards and backwards to find the next SVs
            // int startIndex = startBreakend.getChrPosIndex();
            // endBreakend = findDeletionBreakend(breakendList, startIndex);

            dgData.addSvBreakend(startBreakend, SV_DRIVER_TYPE_DEL);
        }
        else
        {
            endBreakend = minBreakend;
        }

        if(endBreakend != null)
            dgData.addSvBreakend(endBreakend, SV_DRIVER_TYPE_DEL);

        if(startBreakend != null && endBreakend != null)
            dgData.setFullyMatched(true);

        // find straddling LOH for this DEL
        SvBreakend lohStartBreakend = null;
        SvBreakend lohEndBreakend = null;

        LohEvent matchedLohEvent = null;
        if(startBreakend != null && endBreakend != null)
        {
            for (final LohEvent lohEvent : mLohEventList)
            {
                if(!lohEvent.Chromosome.equals(startBreakend.chromosome()))
                    continue;

                // the LOH just needs to straddle one or the other of the min-gene breakends
                if(lohEvent.PosStart < startBreakend.position() && lohEvent.PosEnd > startBreakend.position())
                {
                    SvBreakend lohStart = lohEvent.getBreakend(true);
                    SvBreakend lohEnd = lohEvent.getBreakend(false);

                    if(!lohEvent.isValid() || lohStart == startBreakend || lohEnd == endBreakend)
                    {
                        // the LOH event was likely confused by 2 overlapping DELs
                        if(lohStart == startBreakend && lohEnd != null && lohEnd != endBreakend && lohEnd.getSV().type() == DEL)
                        {
                            lohEndBreakend = lohEnd;
                            lohStartBreakend = lohEnd.getOtherBreakend();
                        }
                        else if(lohEnd == endBreakend && lohStart != null && lohStart != startBreakend && lohStart.getSV().type() == DEL)
                        {
                            lohEndBreakend = lohStart.getOtherBreakend();
                            lohStartBreakend = lohStart;
                        }
                    }
                    else
                    {
                        matchedLohEvent = lohEvent;
                        lohStartBreakend = lohStart;
                        lohEndBreakend = lohEnd;
                    }

                    break;
                }
            }
        }

        // look to next event
        if(lohStartBreakend != null && lohEndBreakend != null)
        {
            dgData.addSvBreakend(lohStartBreakend, SV_DRIVER_TYPE_LOH);
            dgData.addSvBreakend(lohEndBreakend, SV_DRIVER_TYPE_LOH);
        }
        else if(lohStartBreakend != null)
        {
            dgData.addSvBreakend(lohStartBreakend, SV_DRIVER_TYPE_LOH + "_" + matchedLohEvent.SegEnd);
        }
        else if(lohEndBreakend != null)
        {
            dgData.addSvBreakend(lohEndBreakend, SV_DRIVER_TYPE_LOH + "_" + matchedLohEvent.SegStart);
        }
        else if(matchedLohEvent != null)
        {
            dgData.addMatchInfo(String.format("%s;%s", matchedLohEvent.SegStart, matchedLohEvent.SegEnd));
        }
        else
        {
            dgData.setFullyMatched(false);
        }

        writeDriverData(dgData);
                */
    }

}

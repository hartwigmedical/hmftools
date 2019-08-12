package com.hartwig.hmftools.linx.drivers;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.UNKNOWN;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.analysis.SvClassification.isSimpleSingleSV;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.findCentromereBreakendIndex;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatPloidy;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.GAIN;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.GAIN_ARM;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.GAIN_CHR;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_ARM;
import static com.hartwig.hmftools.linx.drivers.DriverEventType.LOH_CHR;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_ARM_SV;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_DEL;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_DM;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_DUP;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_FOLDBACK;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_MAX_PLOIDY;
import static com.hartwig.hmftools.linx.drivers.DriverGeneEvent.SV_DRIVER_TYPE_TI;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.CENTROMERE_CN;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.P_ARM_TELOMERE_CN;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.Q_ARM_TELOMERE_CN;
import static com.hartwig.hmftools.linx.types.ResolvedType.DEL;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS;
import static com.hartwig.hmftools.linx.types.SvCluster.CLUSTER_ANNOT_DM;
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
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.cn.HomLossEvent;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;
import com.hartwig.hmftools.linx.types.ResolvedType;
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
    private List<LohEvent> mLohEventList;
    private Map<String, double[]> mChrCopyNumberMap;
    private Map<String, List<GeneCopyNumber>> mSampleGeneCopyNumberMap; // loaded from file to avoid DB hits on the massive table
    private VisualiserWriter mVisWriter;

    private static final String GCN_DATA_FILE = "gcn_data_file";

    private static final Logger LOGGER = LogManager.getLogger(DriverGeneAnnotator.class);

    public DriverGeneAnnotator(DatabaseAccess dbAccess, SvGeneTranscriptCollection geneTranscriptCollection, final LinxConfig config)
    {
        mDbAccess = dbAccess;
        mGeneTransCache = geneTranscriptCollection;
        mConfig = config;

        mDriverCatalog = Lists.newArrayList();
        mGeneCopyNumberData = Lists.newArrayList();
        mDriverGeneDataList = Lists.newArrayList();
        mDriverOutputList = Lists.newArrayList();
        mSampleGeneCopyNumberMap = new HashMap();
        mChrCopyNumberMap = null;
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

    public void setCopyNumberData(
            final Map<String, double[]> chrCopyNumberMap,
            final List<LohEvent> lohData,
            final List<HomLossEvent> homLossData)
    {
        mLohEventList = lohData;
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

                if (breakendList == null || breakendList.isEmpty())
                {
                    LOGGER.warn("missing breakend list for chromosome({})", dgData.GeneData.Chromosome);
                    writeDriverData(dgData);
                    continue;
                }

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
    }

    private void annotateBiallelicEvent(DriverGeneData dgData)
    {
        // look for an LOH covering any part of the coding region
        if(dgData.TransData.CodingStart == null || dgData.TransData.CodingEnd == null)
            return;

        long codingStart = dgData.TransData.CodingStart;
        long codingEnd = dgData.TransData.CodingEnd;

        for(final LohEvent lohEvent : mLohEventList)
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

        annotateArmAmplification(dgData);
        annotateBreakendAmplification(dgData, breakendList);

        if(dgData.getEvents().isEmpty())
        {
            // otherwise no event
            DriverGeneEvent event = new DriverGeneEvent(GAIN);
            dgData.addEvent(event);
        }
    }

    private void annotateArmAmplification(final DriverGeneData dgData)
    {
        // check for an arm or whole chromosome amplified above the sample ploidy
        double[] cnData = mChrCopyNumberMap.get(dgData.GeneData.Chromosome);

        double telomereCopyNumber = min(cnData[P_ARM_TELOMERE_CN], cnData[Q_ARM_TELOMERE_CN]);

        double armCopyNumber = min(cnData[CENTROMERE_CN],
                dgData.Arm == CHROMOSOME_ARM_P ? cnData[P_ARM_TELOMERE_CN] : cnData[Q_ARM_TELOMERE_CN]);

        double chromosomeCopyNumber = min(cnData[CENTROMERE_CN], telomereCopyNumber);

        if(chromosomeCopyNumber/mSamplePloidy > 2)
        {
            LOGGER.debug("gene({}) AMP gain from chromosomeCopyNumber({}) vs samplePloidy({})",
                    dgData.GeneData.GeneName, formatPloidy(chromosomeCopyNumber), formatPloidy(mSamplePloidy));

            dgData.addEvent(new DriverGeneEvent(GAIN_CHR));
        }
        else if(armCopyNumber/mSamplePloidy > 2)
        {
            LOGGER.debug("gene({}) AMP gain from armCopyNumber({}) vs samplePloidy({})",
                    dgData.GeneData.GeneName, formatPloidy(armCopyNumber), formatPloidy(mSamplePloidy));

            dgData.addEvent(new DriverGeneEvent(GAIN_ARM));
        }
    }

    private void annotateBreakendAmplification(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        // sum up breakend ploidies from telomere to centromere for the gene in question net off ploidy within a cluster
        long transStart = dgData.TransData.TransStart;
        long transEnd = dgData.TransData.TransEnd;

        int centromereIndex = findCentromereBreakendIndex(breakendList, dgData.Arm);
        int startIndex = dgData.Arm == CHROMOSOME_ARM_P ? 0 : centromereIndex;
        int endIndex = dgData.Arm == CHROMOSOME_ARM_P ? centromereIndex : breakendList.size() - 1;

        for(int i = 0; i <= 1; ++i)
        {
            boolean traverseUp = (i == 0);

            Map<SvCluster,List<SvBreakend>> clusterBreakendsMap = Maps.newHashMap();

            int index = traverseUp ? startIndex : endIndex;
            SvBreakend breakend = null;

            while(true)
            {
                if(breakend != null)
                    index += traverseUp ? 1 : -1;

                if(index < startIndex || index > endIndex)
                    break;

                breakend = breakendList.get(index);
                final SvCluster cluster = breakend.getCluster();

                if(cluster.hasLinkingLineElements())
                    continue;

                // make note of any breakend before the gene
                if((traverseUp && breakend.position() > transEnd) || (!traverseUp && breakend.position() < transStart))
                    break;

                if(isSimpleSingleSV(cluster))
                {
                    final SvVarData var = breakend.getSV();
                    if(var.type() != DUP || !traverseUp)
                        continue;

                    if(var.position(true) <= transStart && var.position(false) >= transEnd)
                    {
                        clusterBreakendsMap.put(cluster, Lists.newArrayList(breakend));
                    }

                    continue;
                }

                List<SvBreakend> clusterBreakends = clusterBreakendsMap.get(cluster);

                boolean breakendFacesGene = ((breakend.orientation() != 1) == traverseUp);

                if(clusterBreakends == null)
                {
                    // only start with a breakend facing the gene
                    if(!breakendFacesGene)
                        continue;

                    clusterBreakendsMap.put(cluster, Lists.newArrayList(breakend));
                }
                else
                {
                    clusterBreakends.add(breakend);
                }
            }

            for(Map.Entry<SvCluster,List<SvBreakend>> entry : clusterBreakendsMap.entrySet())
            {
                final List<SvBreakend> clusterBreakends = entry.getValue();

                double netPloidy = clusterBreakends.stream()
                        .mapToDouble(x -> ((x.orientation() != 1) == traverseUp) ? x.ploidy() : -x.ploidy()).sum();

                if(netPloidy < mSamplePloidy)
                    continue;

                final SvCluster cluster = entry.getKey();

                if(dgData.getEvents().stream().anyMatch(x -> x.getCluster() == cluster))
                    continue;

                // report on type of SV which caused the amp
                DriverGeneEvent event = new DriverGeneEvent(GAIN);
                event.setCluster(cluster);

                if(cluster.getSvCount() == 1)
                {
                    final SvVarData var = cluster.getSV(0);
                    final String matchType = cluster.hasAnnotation(CLUSTER_ANNOT_DM) ? SV_DRIVER_TYPE_DM : SV_DRIVER_TYPE_DUP;
                    event.addSvBreakendPair(var.getBreakend(true), var.getBreakend(false), matchType);
                    dgData.addEvent(event);

                    LOGGER.debug("gene({}) cluster({}) net ploidy({}) from DUP({}) type({})",
                            dgData.GeneData.GeneName, cluster.id(), formatPloidy(netPloidy), var.toString(), matchType);

                    continue;
                }

                SvVarData foldback = null;
                SvBreakend maxPloidyBreakend = null;
                SvLinkedPair maxSpanningLink = null;
                double maxSpanningLinkPloidy = 0;

                for(final SvBreakend clusterBreakend : clusterBreakends)
                {
                    boolean breakendFacesGene = ((clusterBreakend.orientation() != 1) == traverseUp);

                    if(!breakendFacesGene)
                        continue;

                    // check for a TI spanning the gene and select the highest ploidy one
                    final SvLinkedPair spanningLink = clusterBreakend.getSV().getLinkedPairs(clusterBreakend.usesStart()).stream()
                            .filter(x -> x.getBreakend(true).position() < transStart && x.getBreakend(false).position() > transEnd)
                            .findFirst().orElse(null);

                    if(spanningLink != null && spanningLink != maxSpanningLink)
                    {
                        // calculate the total ploidy for this link from all chains it is in
                        double linkPloidy = 0;
                        for (final SvChain chain : cluster.getChains())
                        {
                            int linkRepeats = (int) chain.getLinkedPairs().stream().filter(x -> x.matches(spanningLink)).count();
                            linkPloidy += linkRepeats * chain.ploidy();
                        }

                        if(maxSpanningLink == null || maxSpanningLinkPloidy < linkPloidy)
                        {
                            maxSpanningLink = spanningLink;
                            maxSpanningLinkPloidy = linkPloidy;
                        }
                    }

                    if(clusterBreakend.isFoldback() && (foldback == null || clusterBreakend.ploidy() > foldback.ploidy()))
                    {
                        foldback = clusterBreakend.getSV();
                    }

                    if(maxPloidyBreakend == null || clusterBreakend.ploidy() > maxPloidyBreakend.ploidy())
                    {
                        maxPloidyBreakend = clusterBreakend;
                    }
                }

                if(maxSpanningLink != null)
                {
                    LOGGER.debug("gene({}) cluster({}) net ploidy({}) from TI({}) ploidy({})",
                            dgData.GeneData.GeneName, cluster.id(), formatPloidy(netPloidy),
                            maxSpanningLink, formatPloidy(maxSpanningLinkPloidy));

                    event.addSvBreakendPair(maxSpanningLink.getBreakend(true), maxSpanningLink.getBreakend(false), SV_DRIVER_TYPE_TI);
                }
                else if(foldback != null)
                {
                    LOGGER.debug("gene({}) cluster({}) net ploidy({}) from foldback({}) ploidy({})",
                            dgData.GeneData.GeneName, cluster.id(), formatPloidy(netPloidy),
                            foldback, formatPloidy(foldback.ploidy()));

                    event.addSvBreakendPair(foldback.getBreakend(true), foldback.getBreakend(false), SV_DRIVER_TYPE_FOLDBACK);
                }
                else if(maxPloidyBreakend != null)
                {
                    LOGGER.debug("gene({}) cluster({}) net ploidy({}) from max-ploidy breakend({}) ploidy({})",
                            dgData.GeneData.GeneName, cluster.id(), formatPloidy(netPloidy),
                            maxPloidyBreakend, formatPloidy(maxPloidyBreakend.ploidy()));

                    // report highest ploidy facing SV
                    event.addSvBreakendPair(maxPloidyBreakend, null, SV_DRIVER_TYPE_MAX_PLOIDY);
                }

                dgData.addEvent(event);
            }
        }
    }

    private void annotateBreakendAmplification_old(final DriverGeneData dgData, final List<SvBreakend> breakendList)
    {
        // trying to find the breakends which amplify this gene
        // take any DUP or TI which straddles the gene

        // additionally look for any unbalanced breakend ploidy from one or more clusters as any explanation

        long transStart = dgData.TransData.TransStart;
        long transEnd = dgData.TransData.TransEnd;

        SvBreakend foldbackBreakend = null;
        List<Integer> reportedClusters = Lists.newArrayList();

        for (int i = 0; i < breakendList.size(); ++i)
        {
            final SvBreakend breakend = breakendList.get(i);

            if(breakend.position() > transStart)
            {
                // past the gene but AMP may be explained by an unchained foldback
                if(foldbackBreakend == null && breakend.orientation() == 1
                && breakend.getSV().getFoldbackBreakend(breakend.usesStart()) != null)
                {
                    foldbackBreakend = breakend;
                    break;
                }

                continue; // looking for a foldback as an explanation
            }

            if(breakend.orientation() == 1)
                continue;

            final SvVarData varStart = breakend.getSV();

            // only report on each cluster once
            if(reportedClusters.contains(varStart.getCluster().id()))
                continue;

            if(breakend.orientation() == -1 && varStart.getFoldbackBreakend(breakend.usesStart()) != null)
                foldbackBreakend = breakend;

            if(varStart.type() == DUP && varStart.getCluster().getSvCount() == 1)
            {
                if(varStart.position(true) > transStart || varStart.position(false) < transEnd)
                    continue;

                // require the DUP to not be interupted by any other SV within this cluster within this gene region
                SvBreakend upperBreakend = varStart.getBreakend(false);

                String matchType = SV_DRIVER_TYPE_DUP;

                if(varStart.getCluster().hasAnnotation(CLUSTER_ANNOT_DM))
                {
                    LOGGER.debug("gene({}) cluster({}) double minute SV({} {})",
                            dgData.GeneData.GeneName, varStart.getCluster().id(),  varStart.posId(), varStart.type());
                    matchType = SV_DRIVER_TYPE_DM;
                }
                else
                {
                    LOGGER.debug(String.format("gene(%s) cluster(%s) single SV(%s %s) cn(%.2f) cnChg(%.2f)",
                            dgData.GeneData.GeneName, varStart.getCluster().id(), varStart.posId(), varStart.type(),
                            varStart.copyNumber(true), varStart.copyNumberChange(true)));
                }

                DriverGeneEvent event = new DriverGeneEvent(GAIN);
                event.addSvBreakendPair(breakend, upperBreakend, matchType);
                dgData.addEvent(event);

                reportedClusters.add(varStart.getCluster().id());
                continue;
            }

            // skip simple clusters that don't cause amplification
            ResolvedType clusterType = varStart.getCluster().getResolvedType();

            if(clusterType == DEL || clusterType == RECIP_INV || clusterType == RECIP_TRANS)
                continue;

            // look for the first TI which overlaps the gene region
            List<SvLinkedPair> tiPairs = varStart.getLinkedPairs(breakend.usesStart());

            for(SvLinkedPair tiPair : tiPairs)
            {
                if (breakend.position() + tiPair.length() < transEnd)
                    continue;

                final SvVarData varEnd = tiPair.getOtherSV(varStart);
                final SvBreakend beStart = tiPair.getBreakend(true);
                final SvBreakend beEnd = tiPair.getBreakend(false);

                LOGGER.debug(String.format("gene(%s) cluster(%d) SVs start(%s cn=%.2f cnChg=%.2f) end(%s cn=%.2f cnChg=%.2f) in linked pair",
                        dgData.GeneData.GeneName, varStart.getCluster().id(),
                        varStart.posId(), varStart.copyNumber(beStart.usesStart()), varStart.copyNumberChange(beStart.usesStart()),
                        varEnd.posId(), varEnd.copyNumber(beEnd.usesStart()), varEnd.copyNumberChange(beEnd.usesStart())));

                DriverGeneEvent event = new DriverGeneEvent(GAIN);
                event.addSvBreakendPair(beStart, beEnd, SV_DRIVER_TYPE_TI);
                dgData.addEvent(event);

                reportedClusters.add(varStart.getCluster().id());

                break;
            }
        }

        if(!reportedClusters.isEmpty())
            return;

        if(foldbackBreakend != null)
        {
            SvVarData foldbackSv = foldbackBreakend.getSV();
            SvBreakend otherBreakend = foldbackSv.isChainedFoldback() ?
                    foldbackSv.getChainedFoldbackSv().getChainedFoldbackBreakend() : foldbackBreakend.getOtherBreakend();

            if(otherBreakend != null)
            {
                DriverGeneEvent event = new DriverGeneEvent(GAIN);
                event.addSvBreakendPair(foldbackBreakend, otherBreakend, SV_DRIVER_TYPE_FOLDBACK);
                dgData.addEvent(event);
                return;
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
                mFileWriter.write(",SamplePloidy,Chromosome,Arm,MinCN,CentromereCN,TelomereCN");
                mFileWriter.newLine();
            }

            BufferedWriter writer = mFileWriter;

            final DriverCatalog driverGene = dgData.DriverData;
            final EnsemblGeneData geneData = dgData.GeneData;
            double[] cnData = mChrCopyNumberMap.get(geneData.Chromosome);

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
                        : (breakendPair[SE_END] != null ? breakendPair[SE_END].ploidy() : -1);

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

}

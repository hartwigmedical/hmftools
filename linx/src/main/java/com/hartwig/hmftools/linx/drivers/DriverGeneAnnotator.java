package com.hartwig.hmftools.linx.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.PARTIAL_AMP;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.drivers.GeneCopyNumberRegion.calcGeneCopyNumberRegion;
import static com.hartwig.hmftools.common.purple.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.linx.fusion.DisruptionFinder.getDisruptionGeneTranscripts;
import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.DRIVER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.linx.CohortFileInterface;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.TelomereCentromereCnData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.visualiser.file.VisSampleData;

public class DriverGeneAnnotator implements CohortFileInterface
{
    private final EnsemblDataCache mGeneTransCache;
    private final LinxConfig mConfig;

    private final DriverDataCache mDataCache;
    private final AmplificationDrivers mAmpDrivers;
    private final DeletionDrivers mDelDrivers;

    private String mSampleId;
    private final CohortDataWriter mCohortDataWriter;

    private final List<LinxDriver> mDriverOutputList;

    private PerformanceCounter mPerfCounter;

    // references only
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private final VisSampleData mVisSampleData;

    public DriverGeneAnnotator(
            final EnsemblDataCache geneTransCache, final LinxConfig config,
            final CnDataLoader cnDataLoader, final CohortDataWriter cohortDataWriter, final VisSampleData visSampleData)
    {
        mGeneTransCache = geneTransCache;
        mConfig = config;

        mDriverOutputList = Lists.newArrayList();

        mCohortDataWriter = cohortDataWriter;

        mDataCache = new DriverDataCache(cnDataLoader, mGeneTransCache, config.DriverGenes);
        mAmpDrivers = new AmplificationDrivers(mDataCache);

        Map<String,List<String>> disruptionGeneTranscripts = getDisruptionGeneTranscripts(config.DriverGenes, true, geneTransCache);
        mDelDrivers = new DeletionDrivers(disruptionGeneTranscripts, mDataCache);

        mVisSampleData = visSampleData;

        mPerfCounter = new PerformanceCounter("Drivers");
    }

    public PerformanceCounter getPerfCounter() { return mPerfCounter; }

    public void setSamplePurityData(double ploidy, boolean isMale)
    {
        mDataCache.setSamplePurityData(ploidy, isMale);
    }

    public List<DriverGeneData> getDriverGeneDataList() { return mDataCache.getDriverGeneDataList(); }
    public List<DriverCatalog> getDriverCatalog() { return mDataCache.getDriverCatalog(); }
    public List<LinxDriver> getDriverOutputList() { return mDriverOutputList; }

    public void annotateSVs(final String sampleId, final Map<String,List<SvBreakend>> chrBreakendMap)
    {
        mPerfCounter.start();

        mSampleId = sampleId;
        mDataCache.setSampleId(sampleId);
        mChrBreakendMap = chrBreakendMap;
        mDriverOutputList.clear();

        if(!mConfig.PurpleDataPath.isEmpty())
        {
            String samplePurpleDir = mConfig.PurpleDataPath.contains("*") ?
                    mConfig.PurpleDataPath.replaceAll("\\*", sampleId) : mConfig.PurpleDataPath;

            mDataCache.loadDataFromFile(samplePurpleDir);
        }

        // types handled:
        // - TSG biallelic point mutations
        // - AMPs
        // - DELs

        for(final DriverCatalog driverGene : mDataCache.getDriverCatalog())
        {
            if(!driverTypeHandled(driverGene))
                continue;

            final GeneData geneData = mGeneTransCache.getGeneDataByName(driverGene.gene());

            if(geneData == null)
            {
                LNX_LOGGER.warn("driver gene({}) data not found", driverGene.gene());
                continue;
            }

            final TranscriptData transcriptData = mGeneTransCache.getTranscriptData(geneData.GeneId, driverGene.transcript());

            if(transcriptData == null)
            {
                LNX_LOGGER.warn("driver gene({}) trans data not found", driverGene.gene());
                continue;
            }

            GeneCopyNumber geneCN = mDataCache.findGeneCopyNumber(driverGene.gene());
            GeneCopyNumberRegion geneMinCopyNumber;

            if(geneCN != null)
            {
                geneMinCopyNumber = new GeneCopyNumberRegion(geneCN);
            }
            else
            {
                geneMinCopyNumber = calcGeneCopyNumberRegion(transcriptData, mDataCache.CopyNumberData.getChrCnDataMap().get(geneData.Chromosome));

                if(geneMinCopyNumber == null)
                {
                    LNX_LOGGER.warn("sample({}) gene({}) min copy number data not found for driver",
                            mDataCache.sampleId(), driverGene.gene());
                    continue;
                }
            }

            DriverGeneData dgData = new DriverGeneData(driverGene, geneData, transcriptData, geneMinCopyNumber);
            mDataCache.getDriverGeneDataList().add(dgData);

            if(dgData.DriverData.driver() == DriverType.DEL)
            {
                mDelDrivers.annotateDeleteEvent(dgData);
            }
            else if(driverGene.category() == TSG && driverGene.biallelic())
            {
                mDelDrivers.annotateBiallelicEvent(dgData);
            }
            else if(driverGene.driver() == DriverType.AMP || driverGene.driver() == PARTIAL_AMP)
            {
                final List<SvBreakend> breakendList = mChrBreakendMap.get(dgData.GeneInfo.Chromosome);
                mAmpDrivers.annotateAmplification(dgData, breakendList);
            }

            writeDriverData(dgData);
        }

        final List<DriverGeneData> disDelDrivers = mDelDrivers.findDisruptiveDelDrivers(mChrBreakendMap);
        disDelDrivers.forEach(x -> writeDriverData(x));

        mChrBreakendMap = null;

        mPerfCounter.stop();
    }

    private boolean driverTypeHandled(final DriverCatalog driverGene)
    {
        if(driverGene.driver() == DriverType.DEL || driverGene.driver() == DriverType.AMP || driverGene.driver() == PARTIAL_AMP)
            return true;

        if(driverGene.driver() == DriverType.MUTATION && driverGene.category() == TSG && driverGene.biallelic())
        {
            // need to look for an LOH
            return true;
        }

        return false;
    }

    private static final String COHORT_WRITER_DRIVER = "Driver";

    @Override
    public String fileType() { return COHORT_WRITER_DRIVER; }

    @Override
    public BufferedWriter createWriter(final String outputDir)
    {
        try
        {
            String outputFileName = outputDir + "LNX_DRIVERS.csv";

            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("SampleId,Gene,Category,DriverType,LikelihoodMethod,Likelihood");
            writer.write(",FullyMatched,EventType,ClusterId,ClusterCount,ResolvedType");
            writer.write(",Chromosome,Arm,SamplePloidy,GeneMinCN,CentromereCN,TelomereCN,CNGain");
            writer.write(",SvIdStart,SvIdEnd,SvPosStart,SvPosEnd,SvMatchType");
            writer.newLine();
            return writer;
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("failed to initialise cohort driver file: {}", e.toString());
            return null;
        }
    }

    private void writeDriverData(final DriverGeneData dgData)
    {
        if(!dgData.TransData.IsCanonical) // currently ignored since would duplicate everything
            return;

        // convert to a sample driver record
        if(mVisSampleData != null)
        {
            for(final DriverGeneEvent driverEvent : dgData.getEvents())
            {
                int clusterId = driverEvent.getCluster() != null ? driverEvent.getCluster().id() : -1;

                mDriverOutputList.add(ImmutableLinxDriver.builder()
                        .clusterId(clusterId)
                        .gene(dgData.GeneInfo.GeneName)
                        .eventType(driverEvent.Type)
                        .build());

                mVisSampleData.addGeneExonData(clusterId, dgData.GeneInfo.GeneId, dgData.GeneInfo.GeneName,
                        "", 0, dgData.GeneInfo.Chromosome, DRIVER);
            }
        }

        if(!mConfig.hasMultipleSamples())
            return;

        final DriverCatalog driverGene = dgData.DriverData;
        final GeneData geneData = dgData.GeneInfo;

        final TelomereCentromereCnData tcData = mDataCache.CopyNumberData.getChrTeleCentroData().get(dgData.GeneInfo.Chromosome);

        double centromereCopyNumber = 0;
        double telomereCopyNumber = 0;

        if(tcData == null)
        {
            LNX_LOGGER.warn("chromosome({}) missing centro-telo data", dgData.GeneInfo.Chromosome);
        }
        else
        {
            centromereCopyNumber = dgData.Arm == P_ARM ? tcData.CentromerePArm : tcData.CentromereQArm;
            telomereCopyNumber = dgData.Arm == P_ARM ? tcData.TelomerePArm : tcData.TelomereQArm;
        }

        List<String> outputLines = Lists.newArrayList();

        for(final DriverGeneEvent driverEvent : dgData.getEvents())
        {
            final SvBreakend[] breakendPair = driverEvent.getBreakendPair();

            StringBuilder sb = new StringBuilder();

            sb.append(String.format("%s,%s,%s,%s,%s,%.2f,%s,%s",
                    mSampleId, driverGene.gene(), driverGene.category(), driverGene.driver(),
                    driverGene.likelihoodMethod(), driverGene.driverLikelihood(), dgData.fullyMatched(), driverEvent.Type));

            final SvCluster cluster = driverEvent.getCluster();

            if(cluster != null)
            {
                sb.append(String.format(",%d,%d,%s", cluster.id(), cluster.getSvCount(), cluster.getResolvedType()));
            }
            else
            {
                sb.append(String.format(",-1,0,"));
            }

            // breakend info if present

            // Chromosome,Arm,SamplePloidy,GeneMinCN,CentromereCN,TelomereCN,CNGain,SvIdStart,SvIdEnd,SvPosStart,SvPosEnd,SvMatchType

            sb.append(String.format(",%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f",
                    geneData.Chromosome, dgData.Arm, mDataCache.samplePloidy(), dgData.CopyNumberRegion.MinCopyNumber,
                    centromereCopyNumber, telomereCopyNumber, driverEvent.getCopyNumberGain()));

            int posStart = breakendPair[SE_START] != null ? breakendPair[SE_START].position() : 0;
            int posEnd = breakendPair[SE_END] != null ? breakendPair[SE_END].position() : 0;
            String svIdStart = breakendPair[SE_START] != null ? breakendPair[SE_START].getSV().idStr() : "-1";
            String svIdEnd = breakendPair[SE_END] != null ? breakendPair[SE_END].getSV().idStr() : "-1";

            sb.append(String.format(",%s,%s,%d,%d,%s",
                    svIdStart, svIdEnd, posStart, posEnd, driverEvent.getSvInfo()));

            outputLines.add(sb.toString());
        }

        mCohortDataWriter.write(this, outputLines);
    }

    @VisibleForTesting
    public void addDriverGene(final DriverCatalog driver, final GeneCopyNumber geneCopyNumber)
    {
        mDataCache.clearCache();

        mDataCache.getDriverCatalog().add(driver);
        mDataCache.getGeneCopyNumberData().add(geneCopyNumber);
    }

    @VisibleForTesting
    public void clearResults(boolean clearDrivers)
    {
        mDataCache.getDriverGeneDataList().clear();

        if(clearDrivers)
            mDataCache.getDriverCatalog().clear();
    }
}

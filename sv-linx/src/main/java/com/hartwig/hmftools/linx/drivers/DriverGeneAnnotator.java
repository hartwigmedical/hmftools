package com.hartwig.hmftools.linx.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.GERMLINE;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.PARTIAL_AMP;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.drivers.GeneCopyNumberRegion.calcGeneCopyNumberRegion;
import static com.hartwig.hmftools.linx.fusion.DisruptionFinder.disruptionGeneIds;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.DRIVER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.linx.cn.TelomereCentromereCnData;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.visualiser.file.VisualiserWriter;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DriverGeneAnnotator
{
    private final DatabaseAccess mDbAccess;
    private final EnsemblDataCache mGeneTransCache;
    private final LinxConfig mConfig;

    private final DriverDataCache mDataCache;
    private final AmplificationDrivers mAmpDrivers;
    private final DeletionDrivers mDelDrivers;

    private BufferedWriter mFileWriter;

    private String mSampleId;
    private String mOutputDir;

    private final List<LinxDriver> mDriverOutputList;

    private PerformanceCounter mPerfCounter;

    // references only
    private Map<String, List<SvBreakend>> mChrBreakendMap;
    private VisualiserWriter mVisWriter;

    public static final String LINX_DRIVER_CATALOG = ".linx.driver.catalog.tsv";

    public DriverGeneAnnotator(DatabaseAccess dbAccess, final EnsemblDataCache geneTransCache,
            final LinxConfig config, final CnDataLoader cnDataLoader)
    {
        mDbAccess = dbAccess;
        mGeneTransCache = geneTransCache;
        mConfig = config;

        mDriverOutputList = Lists.newArrayList();

        mFileWriter = null;
        mOutputDir = mConfig.OutputDataPath;

        mGeneTransCache.createGeneNameIdMap();

        mDataCache = new DriverDataCache(dbAccess, cnDataLoader, mGeneTransCache);
        mAmpDrivers = new AmplificationDrivers(mDataCache);
        mDelDrivers = new DeletionDrivers(disruptionGeneIds(config.DriverGenes, geneTransCache), mDataCache);

        mVisWriter = null;

        mPerfCounter = new PerformanceCounter("Drivers");
    }

    public void setSamplePurityData(double ploidy, boolean isMale)
    {
        mDataCache.setSamplePurityData(ploidy, isMale);
    }

    public void setVisWriter(VisualiserWriter writer) { mVisWriter = writer; }
    public final List<DriverGeneData> getDriverGeneDataList() { return mDataCache.getDriverGeneDataList(); }

    public void annotateSVs(final String sampleId, final Map<String, List<SvBreakend>> chrBreakendMap)
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
        else if(mDbAccess != null)
        {
            mDataCache.loadDataFromDatabase();
        }

        // types handled:
        // - TSG biallelic point mutations
        // - AMPs
        // - DELs

        for (final DriverCatalog driverGene : mDataCache.getDriverCatalog())
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

            if(canonicalTrans == null)
            {
                List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);
                LNX_LOGGER.warn("driver gene({}) canonical trans not found", driverGene.gene());
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
                geneMinCopyNumber = calcGeneCopyNumberRegion(canonicalTrans, mDataCache.CopyNumberData.getChrCnDataMap().get(geneData.Chromosome));

                if (geneMinCopyNumber == null)
                {
                    LNX_LOGGER.warn("sample({}) gene({}) min copy number data not found for driver",
                            mDataCache.sampleId(), driverGene.gene());
                    continue;
                }
            }

            DriverGeneData dgData = new DriverGeneData(driverGene, geneData, canonicalTrans, geneMinCopyNumber);
            mDataCache.getDriverGeneDataList().add(dgData);

            if(dgData.DriverData.driver() == DriverType.DEL)
            {
                mDelDrivers.annotateDeleteEvent(dgData);
            }
            else if(driverGene.category() == TSG && driverGene.biallelic())
            {
                mDelDrivers.annotateBiallelicEvent(dgData);
            }
            else if (driverGene.driver() == DriverType.AMP || driverGene.driver() == PARTIAL_AMP)
            {
                final List<SvBreakend> breakendList = mChrBreakendMap.get(dgData.GeneData.Chromosome);
                mAmpDrivers.annotateAmplification(dgData, breakendList);
            }

            writeDriverData(dgData);
        }

        final List<DriverGeneData> disDelDrivers = mDelDrivers.findDisruptiveDelDrivers(mChrBreakendMap);
        disDelDrivers.forEach(x -> writeDriverData(x));

        writeSampleDriverData();

        mChrBreakendMap = null;

        mPerfCounter.stop();
    }

    private boolean driverTypeHandled(final DriverCatalog driverGene)
    {
        if(driverGene.driver() == GERMLINE)
            return false;

        if(driverGene.driver() == DriverType.DEL)
            return true;

        if(driverGene.category() == TSG && driverGene.biallelic())
        {
            // need to look for an LOH
            return true;
        }

        if(driverGene.driver() == DriverType.AMP || driverGene.driver() == PARTIAL_AMP)
            return true;

        return false;
    }

    private void writeSampleDriverData()
    {
        if(mConfig.UploadToDB && mDbAccess != null)
        {
            mDbAccess.writeLinxDriverCatalog(mSampleId, mDataCache.getDriverCatalog());
            mDbAccess.writeSvDrivers(mSampleId, mDriverOutputList);
        }

        if(mConfig.hasMultipleSamples() || mOutputDir.isEmpty())
            return;

        // generate an empty Linx driver file even if no annotations were found
        try
        {
            final String driverCatalogFile = mOutputDir + mSampleId + LINX_DRIVER_CATALOG;
            DriverCatalogFile.write(driverCatalogFile, mDataCache.getDriverCatalog());

            final String driversFile = LinxDriver.generateFilename(mOutputDir, mSampleId);
            LinxDriver.write(driversFile, mDriverOutputList);
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
                        "", 0, dgData.GeneData.Chromosome, DRIVER);
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

            final TelomereCentromereCnData tcData = mDataCache.CopyNumberData.getChrTeleCentroData().get(dgData.GeneData.Chromosome);

            double centromereCopyNumber = 0;
            double telomereCopyNumber = 0;

            if(tcData == null)
            {
                LNX_LOGGER.warn("chromosome({}) missing centro-telo data", dgData.GeneData.Chromosome);
            }
            else
            {
                centromereCopyNumber = dgData.Arm == P_ARM ? tcData.CentromerePArm : tcData.CentromereQArm;
                telomereCopyNumber = dgData.Arm == P_ARM ? tcData.TelomerePArm : tcData.TelomereQArm;
            }

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
                        geneData.Chromosome, dgData.Arm, mDataCache.samplePloidy(), dgData.CopyNumberRegion.MinCopyNumber,
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

    public void close()
    {
        if(LNX_LOGGER.isDebugEnabled() || mConfig.hasMultipleSamples())
        {
            mPerfCounter.logStats();
        }

        closeBufferedWriter(mFileWriter);
    }
}

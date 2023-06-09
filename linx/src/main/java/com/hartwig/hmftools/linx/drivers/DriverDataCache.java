package com.hartwig.hmftools.linx.drivers;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_PURPLE_SOMATIC;
import static com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod.DEL;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.cn.CnDataLoader.isMaleSample;
import static com.hartwig.hmftools.linx.drivers.GeneCopyNumberRegion.calcGeneCopyNumberRegion;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.linx.cn.CnDataLoader;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DriverDataCache
{
    public final EnsemblDataCache GeneTransCache;
    public final CnDataLoader CopyNumberData;

    private final DatabaseAccess mDbAccess;

    private final List<DriverCatalog> mDriverCatalog;
    private final List<DriverGeneData> mDriverGeneDataList;

    private final List<GeneCopyNumber> mGeneCopyNumberData;

    private String mSampleId;
    private double mSamplePloidy;
    private boolean mIsMale;

    public DriverDataCache(final DatabaseAccess dbAccess, final CnDataLoader cnDataLoader, final EnsemblDataCache mGeneTransCache)
    {
        CopyNumberData = cnDataLoader;
        GeneTransCache = mGeneTransCache;

        mDbAccess = dbAccess;
        mDriverCatalog = Lists.newArrayList();
        mDriverGeneDataList = Lists.newArrayList();
        mGeneCopyNumberData = Lists.newArrayList();

        mSampleId = "";
        mSamplePloidy = 0;
        mIsMale = false;
    }

    public final List<DriverCatalog> getDriverCatalog() { return mDriverCatalog; }
    public final List<DriverGeneData> getDriverGeneDataList() { return mDriverGeneDataList; }
    public final List<GeneCopyNumber> getGeneCopyNumberData() { return mGeneCopyNumberData; }

    public String sampleId() { return mSampleId; }
    public double samplePloidy() { return mSamplePloidy; }
    public boolean isMale() { return mIsMale; }

    public void setSampleId(final String sampleId) { mSampleId = sampleId; }

    public void clearCache()
    {
        mDriverCatalog.clear();
        mDriverGeneDataList.clear();
        mGeneCopyNumberData.clear();
    }

    public void loadDataFromDatabase()
    {
        clearCache();

        final PurityContext purityContext = mDbAccess.readPurityContext(mSampleId);

        if(purityContext == null)
        {
            LNX_LOGGER.error("sample({}) purity record not found", mSampleId);
            return;
        }

        setSamplePurityData(purityContext.bestFit().ploidy(), isMaleSample(purityContext));

        // add records but filter out any previously added by Linx or any other germline
        mDriverCatalog.addAll(
                mDbAccess.readDriverCatalog(mSampleId).stream()
                        .filter(x -> DRIVERS_PURPLE_SOMATIC.contains(x.driver()))
                        .collect(Collectors.toList()));

        LNX_LOGGER.debug("retrieved {} driver gene records", mDriverCatalog.size());

        if(!mDriverCatalog.isEmpty())
        {
            final List<String> driverGenes = mDriverCatalog.stream().map(x -> x.gene()).collect(Collectors.toList());
            mGeneCopyNumberData.addAll(mDbAccess.readGeneCopynumbers(mSampleId, driverGenes));
        }
    }

    public void loadDataFromFile(final String purpleDataPath)
    {
        clearCache();

        try
        {
            final PurityContext purityContext = PurityContextFile.read(purpleDataPath, mSampleId);
            setSamplePurityData(purityContext.bestFit().ploidy(), isMaleSample(purityContext));

            mDriverCatalog.addAll(
                    DriverCatalogFile.read(DriverCatalogFile.generateSomaticFilename(purpleDataPath, mSampleId)).stream()
                            .filter(x -> DRIVERS_PURPLE_SOMATIC.contains(x.driver()))
                            .collect(Collectors.toList()));

            mGeneCopyNumberData.addAll(
                    GeneCopyNumberFile.read(GeneCopyNumberFile.generateFilenameForReading(purpleDataPath, mSampleId)));
        }
        catch(IOException e)
        {
            LNX_LOGGER.error("failed to load driver catalog or purity context: {}", e.toString());
        }
    }

    public void setSamplePurityData(double ploidy, boolean isMale)
    {
        mSamplePloidy = ploidy;
        mIsMale = isMale;
    }

    public DriverGeneData createDriverData(final BreakendGeneData gene, final TranscriptData transData, final DriverType driverType)
    {
        GeneCopyNumber gcnData = null;

        if(mDbAccess != null)
        {
            final List<GeneCopyNumber> geneCopyNumbers = mDbAccess.readGeneCopynumbers(mSampleId, Lists.newArrayList(gene.geneName()));

            if(geneCopyNumbers != null)
            {
                for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
                {
                    if(geneCopyNumber.transName().equals(transData.TransName))
                    {
                        gcnData = geneCopyNumber;
                        break;
                    }
                }

                if(gcnData == null && !geneCopyNumbers.isEmpty())
                    gcnData = geneCopyNumbers.get(0);
            }
        }
        else
        {
            gcnData = findGeneCopyNumber(gene.geneName());
        }

        final DriverCatalog driverRecord = ImmutableDriverCatalog.builder()
                .driver(driverType)
                .category(TSG)
                .gene(gene.geneName())
                .transcript(transData.TransName)
                .isCanonical(transData.IsCanonical)
                .chromosome(gene.chromosome())
                .chromosomeBand(gene.GeneData.KaryotypeBand)
                .likelihoodMethod(DEL)
                .driverLikelihood(1.0)
                .missense(0)
                .nonsense(0)
                .splice(0)
                .inframe(0)
                .frameshift(0)
                .biallelic(true)
                .minCopyNumber(gcnData != null ? gcnData.minCopyNumber() : 0)
                .maxCopyNumber(gcnData != null ? gcnData.maxCopyNumber() : 0)
                .build();

        mDriverCatalog.add(driverRecord);

        final GeneData geneData = GeneTransCache.getGeneDataByName(gene.geneName());

        if(geneData == null)
        {
            LNX_LOGGER.warn("gene({}) no Ensembl gene data found", gene.geneName());
            return null;
        }

        GeneCopyNumberRegion copyNumberRegion;

        if(gcnData != null)
            copyNumberRegion = new GeneCopyNumberRegion(gcnData);
        else
            copyNumberRegion = calcGeneCopyNumberRegion(transData, CopyNumberData.getChrCnDataMap().get(geneData.Chromosome));

        DriverGeneData dgData = new DriverGeneData(driverRecord, geneData, transData, copyNumberRegion);
        mDriverGeneDataList.add(dgData);

        return dgData;
    }

    public final GeneCopyNumber findGeneCopyNumber(final String geneName)
    {
        for(GeneCopyNumber geneCN : mGeneCopyNumberData)
        {
            if(geneCN.geneName().equals(geneName))
                return geneCN;
        }

        return null;
    }

}

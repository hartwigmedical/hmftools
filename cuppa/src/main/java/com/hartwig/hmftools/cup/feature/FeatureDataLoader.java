package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus.MSS;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_CHROMOSOME;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_AMP;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_DEL;
import static com.hartwig.hmftools.cup.feature.ViralInsertionType.OTHER;
import static com.hartwig.hmftools.cup.feature.ViralInsertionType.fromVirusName;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRIVERCATALOG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVFUSION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.VIRALINSERTION;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContextFile;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.variant.structural.linx.FusionPhasedType;
import com.hartwig.hmftools.common.variant.structural.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.cup.somatics.SomaticDataLoader;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.jooq.Record;
import org.jooq.Result;

public class FeatureDataLoader
{
    public static boolean loadDriversFromCohortFile(final String filename, final Map<String,List<SampleFeatureData>> sampleDrivers)
    {
        if(filename == null)
            return true;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());
            fileData.remove(0);

            for(final String line : fileData)
            {
                final SampleFeatureData driverData = SampleFeatureData.from(line);

                if(driverData != null)
                {
                    List<SampleFeatureData> drivers = sampleDrivers.get(driverData.SampleId);

                    if(drivers == null)
                    {
                        sampleDrivers.put(driverData.SampleId, Lists.newArrayList(driverData));
                    }
                    else
                    {
                        drivers.add(driverData);
                    }
                }
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read sample driver data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public static boolean loadFeaturesFromFile(
            final String sampleId, final String sampleDataDir, final String sampleVcfFile,
            final Map<String,List<SampleFeatureData>> sampleFeaturesMap)
    {
        try
        {
            String viralInsertFilename = LinxViralInsertion.generateFilename(sampleDataDir, sampleId);

            final List<LinxViralInsertion> viralInserts = Files.exists(Paths.get(viralInsertFilename)) ?
                    LinxViralInsertion.read(viralInsertFilename) : Lists.newArrayList();

            final String fusionsFilename = LinxFusion.generateFilename(sampleDataDir, sampleId);

            final List<LinxFusion> fusions = Files.exists(Paths.get(fusionsFilename)) ?
                    LinxFusion.read(fusionsFilename) : Lists.newArrayList();

            final String purpleDriverCatalogFilename = DriverCatalogFile.generateSomaticFilenameForReading(sampleDataDir, sampleId);

            final List<DriverCatalog> drivers = Files.exists(Paths.get(purpleDriverCatalogFilename)) ?
                    DriverCatalogFile.read(purpleDriverCatalogFilename) : Lists.newArrayList();

            final List<String> indelGenes = checkIndels(sampleId, sampleDataDir, null) ?
                    loadSpecificIndels(sampleId, sampleVcfFile) : null;

            mapFeatureData(sampleId, sampleFeaturesMap, drivers, fusions, viralInserts, indelGenes);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to load drivers, fusion and virus data files: {}", e.toString());
            return false;
        }

        return true;
    }

    private static boolean checkIndels(final String sampleId, final String sampleDataDir, final DatabaseAccess dbAccess)
    {
        PurityContext purityContext = null;

        if(dbAccess != null)
        {
            purityContext = dbAccess.readPurityContext(sampleId);
        }
        else
        {
            try
            {
                purityContext = PurityContextFile.read(sampleDataDir, sampleId);
            }
            catch (Exception e)
            {
                CUP_LOGGER.error("sample({}) check indels - failed to load purity file( from dir{}): {}",
                        sampleId, sampleDataDir, e.toString());
                return false;
            }
        }

        return purityContext != null && purityContext.microsatelliteStatus() == MSS;
    }

    public static boolean loadFeaturesFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds,
            final Map<String,List<SampleFeatureData>> sampleFeaturesMap, boolean skipZeroDriverLikelihood)
    {
        if(dbAccess == null)
            return false;

        final String specificSampleId = sampleIds.size() == 1 ? sampleIds.get(0) : null;

        final Map<String,List<DriverCatalog>> sampleDriverMap = getAllDrivers(dbAccess, skipZeroDriverLikelihood, specificSampleId);

        final Map<String,List<LinxFusion>> sampleFusionMap = getAllFusions(dbAccess, specificSampleId);

        final Map<String,List<LinxViralInsertion>> sampleVirusMap = getAllViruses(dbAccess, specificSampleId);

        final Map<String,List<String>> sampleIndelMap = getAllIndels(dbAccess, specificSampleId);

        int i = 0;
        int nextLog = 100;

        for(final String sampleId : sampleIds)
        {
            final List<DriverCatalog> drivers = sampleDriverMap.get(sampleId);
            final List<LinxFusion> fusions = sampleFusionMap.get(sampleId);
            final List<LinxViralInsertion> viralInserts = sampleVirusMap.get(sampleId);

            final List<String> indelGenes = checkIndels(sampleId, null, dbAccess) ? sampleIndelMap.get(sampleId) : null;

            mapFeatureData(sampleId, sampleFeaturesMap, drivers, fusions, viralInserts, indelGenes);

            ++i;
            if(i >= nextLog)
            {
                nextLog += 100;
                CUP_LOGGER.debug("loaded {} sample feature data sets", i);
            }
        }

        return true;
    }

    private static final Map<String,List<LinxFusion>> getAllFusions(final DatabaseAccess dbAccess, final String specificSampleId)
    {
        final Map<String,List<LinxFusion>> sampleFusionMap = Maps.newHashMap();

        Result<Record> result = dbAccess.context().select().from(SVFUSION)
                .where(SVFUSION.REPORTED.eq((byte)1))
                .and(specificSampleId != null ? SVFUSION.SAMPLEID.eq(specificSampleId) : SVFUSION.SAMPLEID.isNotNull())
                .fetch();

        for (Record record : result)
        {
            final String sampleId = record.getValue(SVFUSION.SAMPLEID);

            LinxFusion fusion = ImmutableLinxFusion.builder()
                    .fivePrimeBreakendId(record.getValue(SVFUSION.FIVEPRIMEBREAKENDID).intValue())
                    .threePrimeBreakendId(record.getValue(SVFUSION.THREEPRIMEBREAKENDID).intValue())
                    .name(record.getValue(SVFUSION.NAME))
                    .reported(record.getValue(SVFUSION.REPORTED) == 1)
                    .reportedType(record.getValue(SVFUSION.REPORTEDTYPE))
                    .likelihood(FusionLikelihoodType.valueOf(record.getValue(SVFUSION.LIKELIHOOD)))
                    .phased(FusionPhasedType.valueOf(record.getValue(SVFUSION.PHASED)))
                    .chainLength(record.getValue(SVFUSION.CHAINLENGTH))
                    .chainLinks(record.getValue(SVFUSION.CHAINLINKS))
                    .chainTerminated(record.getValue(SVFUSION.CHAINTERMINATED) == 1)
                    .domainsKept(record.getValue(SVFUSION.DOMAINSKEPT))
                    .domainsLost(record.getValue(SVFUSION.DOMAINSLOST))
                    .skippedExonsUp(record.getValue(SVFUSION.SKIPPEDEXONSUP))
                    .skippedExonsDown(record.getValue(SVFUSION.SKIPPEDEXONSDOWN))
                    .fusedExonUp(record.getValue(SVFUSION.FUSEDEXONUP))
                    .fusedExonDown(record.getValue(SVFUSION.FUSEDEXONDOWN))
                    .geneStart("")
                    .geneContextStart("")
                    .geneTranscriptStart("")
                    .geneEnd("")
                    .geneContextEnd("")
                    .geneTranscriptEnd("")
                    .junctionCopyNumber(0.0)
                    .build();

            final List<LinxFusion> fusions = sampleFusionMap.get(sampleId);
            if(fusions == null)
                sampleFusionMap.put(sampleId, Lists.newArrayList(fusion));
            else
                fusions.add(fusion);
        }

        return sampleFusionMap;
    }

    private static final Map<String,List<LinxViralInsertion>> getAllViruses(final DatabaseAccess dbAccess, final String specificSampleId)
    {
        final Map<String,List<LinxViralInsertion>> sampleVirusMap = Maps.newHashMap();

        Result<Record> result = dbAccess.context().select()
                .from(VIRALINSERTION)
                .where(specificSampleId != null ? VIRALINSERTION.SAMPLEID.eq(specificSampleId) : VIRALINSERTION.SAMPLEID.isNotNull())
                .fetch();

        for (Record record : result)
        {
            final String sampleId = record.getValue(VIRALINSERTION.SAMPLEID);

            LinxViralInsertion viralInsertion = new LinxViralInsertion(
                    sampleId, record.getValue(VIRALINSERTION.SVID),
                    record.getValue(VIRALINSERTION.VIRUSID), record.getValue(VIRALINSERTION.VIRUSNAME));

            final List<LinxViralInsertion> viralInsertions = sampleVirusMap.get(sampleId);
            if(viralInsertions == null)
                sampleVirusMap.put(sampleId, Lists.newArrayList(viralInsertion));
            else
                viralInsertions.add(viralInsertion);
        }

        return sampleVirusMap;
    }

    private static final Map<String,List<String>> getAllIndels(final DatabaseAccess dbAccess, final String specificSampleId)
    {
        final Map<String,List<String>> sampleIndelMap = Maps.newHashMap();

        Result<Record> result = dbAccess.context().select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(specificSampleId != null ? SOMATICVARIANT.SAMPLEID.eq(specificSampleId) : SOMATICVARIANT.SAMPLEID.isNotNull())
                .and(SOMATICVARIANT.GENE.in(INDEL_ALB, INDEL_SFTPB, INDEL_SLC34A2))
                .and(SOMATICVARIANT.REPEATCOUNT.lessOrEqual(INDEL_MAX_REPEAT_COUNT))
                .and(SOMATICVARIANT.TYPE.eq(VariantType.INDEL.toString()))
                .fetch();

        for (Record record : result)
        {
            final String sampleId = record.getValue(SOMATICVARIANT.SAMPLEID);
            final String gene = record.getValue(SOMATICVARIANT.GENE);

            final List<String> genes = sampleIndelMap.get(sampleId);
            if(genes == null)
                sampleIndelMap.put(sampleId, Lists.newArrayList(gene));
            else
                genes.add(gene);
        }

        return sampleIndelMap;
    }

    private static final Map<String,List<DriverCatalog>> getAllDrivers(
            final DatabaseAccess dbAccess, boolean skipZeroDriverLikelihood, final String specificSampleId)
    {
        final Map<String,List<DriverCatalog>> sampleDriverMap = Maps.newHashMap();

        final Result<Record> result = dbAccess.context().select().from(DRIVERCATALOG)
                .where(specificSampleId != null ? DRIVERCATALOG.SAMPLEID.eq(specificSampleId) : DRIVERCATALOG.SAMPLEID.isNotNull())
                .fetch();

        for (Record record : result)
        {
            final String sampleId = record.getValue(DRIVERCATALOG.SAMPLEID);

            DriverCatalog driverCatalog = ImmutableDriverCatalog.builder()
                    .gene(record.getValue(DRIVERCATALOG.GENE))
                    .chromosome(record.getValue(DRIVERCATALOG.CHROMOSOME))
                    .chromosomeBand(record.getValue(DRIVERCATALOG.CHROMOSOMEBAND))
                    .driver(DriverType.valueOf(record.getValue(DRIVERCATALOG.DRIVER)))
                    .category(DriverCategory.valueOf(record.getValue(DRIVERCATALOG.CATEGORY)))
                    .likelihoodMethod(LikelihoodMethod.valueOf(record.getValue(DRIVERCATALOG.LIKELIHOODMETHOD)))
                    .driverLikelihood(record.getValue(DRIVERCATALOG.DRIVERLIKELIHOOD))
                    .dndsLikelihood(record.getValue(DRIVERCATALOG.DNDSLIKELIHOOD))
                    .missense(record.getValue(DRIVERCATALOG.MISSENSE))
                    .nonsense(record.getValue(DRIVERCATALOG.NONSENSE))
                    .splice(record.getValue(DRIVERCATALOG.SPLICE))
                    .inframe(record.getValue(DRIVERCATALOG.INFRAME))
                    .frameshift(record.getValue(DRIVERCATALOG.FRAMESHIFT))
                    .biallelic(record.getValue(DRIVERCATALOG.BIALLELIC) != 0)
                    .minCopyNumber(record.getValue(DRIVERCATALOG.MINCOPYNUMBER))
                    .maxCopyNumber(record.getValue(DRIVERCATALOG.MAXCOPYNUMBER))
                    .build();

            if(skipZeroDriverLikelihood && driverCatalog.driverLikelihood() <= 0)
                continue;

            final List<DriverCatalog> drivers = sampleDriverMap.get(sampleId);
            if(drivers == null)
                sampleDriverMap.put(sampleId, Lists.newArrayList(driverCatalog));
            else
                drivers.add(driverCatalog);
        }

        return sampleDriverMap;
    }

    private static final int INDEL_MAX_REPEAT_COUNT = 6;

    private static List<String> loadSpecificIndels(final String sampleId, final String vcfFile)
    {
        final List<SomaticVariant> variants = SomaticDataLoader.loadSomaticVariants(sampleId, vcfFile, Lists.newArrayList(VariantType.INDEL));

        return variants.stream()
                .filter(x -> x.repeatCount() <= INDEL_MAX_REPEAT_COUNT)
                .filter(x -> x.filter().equals("PASS"))
                .map(x -> x.gene())
                .collect(Collectors.toList());
    }

    private static void mapFeatureData(
            final String sampleId, final Map<String,List<SampleFeatureData>> sampleDrivers, final List<DriverCatalog> drivers,
            final List<LinxFusion> fusions, final List<LinxViralInsertion> viralInserts, final List<String> indelGenes)
    {
        final List<SampleFeatureData> featuresList = Lists.newArrayList();

        if(drivers != null)
        {
            for(final DriverCatalog driver : drivers)
            {
                SampleFeatureData feature = new SampleFeatureData(sampleId, driver.gene(), DRIVER, driver.driverLikelihood());

                if(driver.driver() == DriverType.AMP)
                    feature.ExtraInfo.put(DRIVER_TYPE, DRIVER_TYPE_AMP);
                else if(driver.driver() == DriverType.DEL)
                    feature.ExtraInfo.put(DRIVER_TYPE, DRIVER_TYPE_DEL);

                feature.ExtraInfo.put(DRIVER_CHROMOSOME, driver.chromosome());

                featuresList.add(feature);
            }
        }

        if(fusions != null)
        {
            for(final LinxFusion fusion : fusions)
            {
                if(!fusion.reported())
                    continue;

                if(fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_5.toString())
                || fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_3.toString()))
                {
                    if(fusion.likelihood() != FusionLikelihoodType.HIGH)
                        continue;

                    final String[] genes = fusion.name().split("_");

                    final String fusionName = fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_5.toString()) ?
                            genes[0] + "_PROM5" : genes[1] + "_PROM3";

                    featuresList.add(new SampleFeatureData(sampleId, fusionName, FeatureType.FUSION, 1));
                }
                else
                {
                    featuresList.add(new SampleFeatureData(sampleId, fusion.name(), FeatureType.FUSION, 1));
                }
            }
        }

        if(viralInserts != null)
        {
            final List<SampleFeatureData> viralInsertDataList = Lists.newArrayList();

            viralInserts.stream()
                    .map(x -> new SampleFeatureData(sampleId, fromVirusName(x.VirusName).toString(), FeatureType.VIRUS, 1))
                    .filter(x -> !x.Name.equals(OTHER.toString()))
                    .filter(x -> viralInsertDataList.stream().noneMatch(y -> y.Name.equals(x.Name))) // check for duplicates
                    .forEach(x -> viralInsertDataList.add(x));

            featuresList.addAll(viralInsertDataList);
        }

        if(indelGenes != null)
        {
            final List<SampleFeatureData> indelFeatures = indelGenes.stream()
                    .map(x -> new SampleFeatureData(sampleId, String.format("INDEL_%s", x), FeatureType.INDEL, 1))
                    .filter(x -> !x.Name.equals(OTHER.toString()))
                    .collect(Collectors.toList());

            featuresList.addAll(indelFeatures);
        }

        sampleDrivers.put(sampleId, featuresList);
    }

    private static final String INDEL_ALB = "ALB";
    private static final String INDEL_SFTPB = "SFTPB";
    private static final String INDEL_SLC34A2 = "SLC34A2";

    public static boolean loadRefPrevalenceData(
            final String filename, final Map<String,FeaturePrevCounts> genePrevalenceTotals,
            final Map<String,List<FeaturePrevData>> cancerDriverPrevalence)
    {
        if(filename == null || filename.isEmpty())
            return false;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            fileData.remove(0);

            for(final String line : fileData)
            {
                final FeaturePrevData prevData = FeaturePrevData.from(line);

                if(prevData == null)
                    continue;

                FeaturePrevCounts genePrevTotals = genePrevalenceTotals.get(prevData.Name);

                if(genePrevTotals == null)
                {
                    genePrevTotals = new FeaturePrevCounts();
                    genePrevalenceTotals.put(prevData.Name, genePrevTotals);
                }

                genePrevTotals.MaxPrevalence = max(genePrevTotals.MaxPrevalence, prevData.Prevalence);

                final List<FeaturePrevData> dataList = cancerDriverPrevalence.get(prevData.CancerType);
                if(dataList == null)
                {
                    cancerDriverPrevalence.put(prevData.CancerType, Lists.newArrayList(prevData));
                }
                else
                {
                    dataList.add(prevData);
                }
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read driver prevalence data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

    public static boolean loadRefCancerFeatureAvg(final String filename, final Map<String,Double> cancerFeatureAvgs)
    {
        if(filename == null || filename.isEmpty())
            return false;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            fileData.remove(0);

            for(final String line : fileData)
            {
                final String[] items = line.split(DATA_DELIM);
                final String cancerType = items[0];
                double average = Double.parseDouble(items[1]);
                cancerFeatureAvgs.put(cancerType, average);
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read feature averages data file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }

}

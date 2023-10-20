package com.hartwig.hmftools.cup.feature;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus.MSS;
import static com.hartwig.hmftools.common.virus.VirusLikelihoodType.HIGH;
import static com.hartwig.hmftools.common.virus.VirusLikelihoodType.UNKNOWN;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.CuppaConfig.DATA_DELIM;
import static com.hartwig.hmftools.cup.CuppaRefFiles.purpleSomaticVcfFile;
import static com.hartwig.hmftools.cup.common.CupConstants.KNOWN_MUTATIONS;
import static com.hartwig.hmftools.cup.feature.FeatureType.AMP;
import static com.hartwig.hmftools.cup.feature.FeatureType.DRIVER;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.AMP_CN;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_CHROMOSOME;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_AMP;
import static com.hartwig.hmftools.cup.feature.SampleFeatureData.DRIVER_TYPE_DEL;
import static com.hartwig.hmftools.cup.feature.ViralInsertionType.OTHER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRIVERCATALOG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SVFUSION;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.VIRUSANNOTATION;

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
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.linx.FusionPhasedType;
import com.hartwig.hmftools.common.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.common.virus.ImmutableAnnotatedVirus;
import com.hartwig.hmftools.common.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.common.virus.VirusType;
import com.hartwig.hmftools.common.virus.VirusLikelihoodType;
import com.hartwig.hmftools.cup.somatics.SomaticDataLoader;
import com.hartwig.hmftools.cup.somatics.SomaticVariant;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

import org.apache.logging.log4j.util.Strings;
import org.jooq.Record;
import org.jooq.Result;

public class FeatureDataLoader
{
    public static boolean loadFeaturesFromCohortFile(final String filename, final Map<String,List<SampleFeatureData>> sampleDrivers)
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
            final String sampleId, final String linxDataDir, final String purpleDataDir, final String virusDataDir,
            final Map<String,List<SampleFeatureData>> sampleFeaturesMap)
    {
        // extract features from standard pipeline output files and fail if any cannot be loaded
        try
        {
            String viralAnnotationFilename = AnnotatedVirusFile.generateFileName(virusDataDir, sampleId);

            final List<AnnotatedVirus> virusAnnotations = Lists.newArrayList();

            if(Files.exists(Paths.get(viralAnnotationFilename)))
            {
                AnnotatedVirusFile.read(viralAnnotationFilename).stream().filter(x -> x.reported()).forEach(x -> virusAnnotations.add(x));
            }
            else
            {
                CUP_LOGGER.warn("sample({}) missing viral annotations file: {}", sampleId, viralAnnotationFilename);
            }

            final String fusionsFilename = LinxFusion.generateFilename(linxDataDir, sampleId);

            final List<LinxFusion> fusions = LinxFusion.read(fusionsFilename);

            // load linx drivers if available, otherwise the purple somatic drivers
            final List<DriverCatalog> drivers = Lists.newArrayList();

            final String linxDriverCatalogFilename = LinxDriver.generateCatalogFilename(linxDataDir, sampleId, true);
            final String purpleDriverCatalogFilename = DriverCatalogFile.generateSomaticFilename(purpleDataDir, sampleId);

            if(Files.exists(Paths.get(linxDriverCatalogFilename)))
            {
                drivers.addAll(DriverCatalogFile.read(linxDriverCatalogFilename));
            }
            else if(Files.exists(Paths.get(purpleDriverCatalogFilename)))
            {
                drivers.addAll(DriverCatalogFile.read(purpleDriverCatalogFilename));
            }
            else
            {
                CUP_LOGGER.error("sample({}) failed to load drivers", sampleId);
                return false;
            }

            boolean checkIndels = checkIndels(sampleId, purpleDataDir, null);
            final String somaticVcf = purpleSomaticVcfFile(purpleDataDir, sampleId);
            final List<String> indelGenes = loadSpecificMutations(somaticVcf, checkIndels);

            mapFeatureData(sampleId, sampleFeaturesMap, drivers, fusions, virusAnnotations, indelGenes);
        }
        catch(IOException e)
        {
            CUP_LOGGER.error("sample({}) failed to load drivers, fusion and virus data files: {}", sampleId, e.toString());
            return false;
        }

        return true;
    }

    private static boolean checkIndels(final String sampleId, final String purpleDataDir, final DatabaseAccess dbAccess)
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
                purityContext = PurityContextFile.read(purpleDataDir, sampleId);
            }
            catch (Exception e)
            {
                CUP_LOGGER.error("sample({}) check indels - failed to load purity file( from dir{}): {}",
                        sampleId, purpleDataDir, e.toString());
                return false;
            }
        }

        return purityContext != null && purityContext.microsatelliteStatus() == MSS;
    }

    public static boolean loadFeaturesFromDatabase(
            final DatabaseAccess dbAccess, final List<String> sampleIds, final Map<String,List<SampleFeatureData>> sampleFeaturesMap)
    {
        if(dbAccess == null)
            return false;

        final String specificSampleId = sampleIds.size() == 1 ? sampleIds.get(0) : null;

        final Map<String,List<DriverCatalog>> sampleDriverMap = getAllDrivers(dbAccess, specificSampleId);

        final Map<String,List<LinxFusion>> sampleFusionMap = getAllFusions(dbAccess, specificSampleId);

        final Map<String,List<AnnotatedVirus>> sampleVirusMap = getAllViruses(dbAccess, specificSampleId);

        final Map<String,List<String>> sampleIndelMap = getSpecificMutations(dbAccess, specificSampleId, true);

        int i = 0;
        int nextLog = 100;

        for(final String sampleId : sampleIds)
        {
            final List<DriverCatalog> drivers = sampleDriverMap.get(sampleId);
            final List<LinxFusion> fusions = sampleFusionMap.get(sampleId);
            final List<AnnotatedVirus> virusAnnotations = sampleVirusMap.get(sampleId);

            final List<String> mutationGenes = sampleIndelMap.get(sampleId);

            if(mutationGenes != null && !checkIndels(sampleId, null, dbAccess))
            {
                final List<String> nonIndelMutations = mutationGenes.stream()
                        .filter(x -> !isKnownIndelGene(x))
                        .collect(Collectors.toList());

                mutationGenes.clear();
                mutationGenes.addAll(nonIndelMutations);
            }

            mapFeatureData(sampleId, sampleFeaturesMap, drivers, fusions, virusAnnotations, mutationGenes);

            ++i;
            if(i >= nextLog)
            {
                nextLog += 100;
                CUP_LOGGER.debug("loaded {} sample feature data sets", i);
            }
        }

        return true;
    }

    private static Map<String,List<LinxFusion>> getAllFusions(final DatabaseAccess dbAccess, final String specificSampleId)
    {
        final Map<String,List<LinxFusion>> sampleFusionMap = Maps.newHashMap();

        Result<Record> result = dbAccess.context().select().from(SVFUSION)
                .where(SVFUSION.REPORTED.eq((byte)1))
                .and(specificSampleId != null ? SVFUSION.SAMPLEID.eq(specificSampleId) : SVFUSION.SAMPLEID.isNotNull())
                .fetch();

        for(Record record : result)
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
                    .reportableReasons("")
                    .build();

            final List<LinxFusion> fusions = sampleFusionMap.get(sampleId);
            if(fusions == null)
                sampleFusionMap.put(sampleId, Lists.newArrayList(fusion));
            else
                fusions.add(fusion);
        }

        return sampleFusionMap;
    }

    private static Map<String,List<AnnotatedVirus>> getAllViruses(final DatabaseAccess dbAccess, final String specificSampleId)
    {
        final Map<String,List<AnnotatedVirus>> sampleVirusMap = Maps.newHashMap();

        Result<Record> result = dbAccess.context().select()
                .from(VIRUSANNOTATION)
                .where(specificSampleId != null ? VIRUSANNOTATION.SAMPLEID.eq(specificSampleId) : VIRUSANNOTATION.SAMPLEID.isNotNull())
                .and(VIRUSANNOTATION.REPORTED.eq((byte)1))
                .fetch();

        for (Record record : result)
        {
            final String sampleId = record.getValue(VIRUSANNOTATION.SAMPLEID);

            String interpretation = record.getValue(VIRUSANNOTATION.INTERPRETATION);

            AnnotatedVirus annotatedVirus = ImmutableAnnotatedVirus.builder()
                        .taxid(record.getValue(VIRUSANNOTATION.TAXID))
                        .name(record.getValue(VIRUSANNOTATION.VIRUSNAME))
                        .qcStatus(VirusBreakendQCStatus.valueOf(record.getValue(VIRUSANNOTATION.QCSTATUS)))
                        .interpretation(interpretation == null ? null : VirusType.fromVirusName(interpretation))
                        .reported(record.getValue(VIRUSANNOTATION.REPORTED) == 1)
                        .integrations(record.getValue(VIRUSANNOTATION.INTEGRATIONS))
                        .percentageCovered(record.getValue(VIRUSANNOTATION.PERCENTAGECOVERED))
                        .meanCoverage(record.getValue(VIRUSANNOTATION.MEANCOVERAGE))
                        .expectedClonalCoverage(record.getValue(VIRUSANNOTATION.EXPECTEDCLONALCOVERAGE))
                        .percentageCovered(record.getValue(VIRUSANNOTATION.PERCENTAGECOVERED))
                        .virusDriverLikelihoodType(VirusLikelihoodType.valueOf(record.getValue(VIRUSANNOTATION.LIKELIHOOD)))
                        .build();

            List<AnnotatedVirus> annotatedVirusList = sampleVirusMap.get(sampleId);
            if(annotatedVirusList == null)
            {
                annotatedVirusList = Lists.newArrayList();
                sampleVirusMap.put(sampleId, annotatedVirusList);
            }

            annotatedVirusList.add(annotatedVirus);
        }

        return sampleVirusMap;
    }

    private static final int INDEL_MAX_REPEAT_COUNT = 6;
    private static final String INDEL_ALB = "ALB";
    private static final String INDEL_SFTPB = "SFTPB";
    private static final String INDEL_SLC34A2 = "SLC34A2";

    private static boolean isKnownIndelGene(final String gene)
    {
        return gene.equals(INDEL_ALB) || gene.equals(INDEL_SFTPB) || gene.equals(INDEL_SLC34A2);
    }

    public static boolean isKnownIndel(final String gene, final int repeatCount, final VariantType variantType)
    {
        return isKnownIndelGene(gene) && variantType == INDEL && repeatCount <= INDEL_MAX_REPEAT_COUNT;
    }

    private static final String MUTATION_EGFR = "EGFR";

    private static boolean isKnownMutation(
            final String gene, final VariantType variantType, final int position, final String ref, final String alt)
    {
        return KNOWN_MUTATIONS.stream().anyMatch(x -> x.matches(gene, variantType, ref, alt, position));
    }

    private static Map<String,List<String>> getSpecificMutations(
            final DatabaseAccess dbAccess, final String specificSampleId, boolean checkIndels)
    {
        final Map<String,List<String>> sampleMutationMap = Maps.newHashMap();

        Result<Record> result = dbAccess.context().select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(specificSampleId != null ? SOMATICVARIANT.SAMPLEID.eq(specificSampleId) : SOMATICVARIANT.SAMPLEID.isNotNull())
                .and(SOMATICVARIANT.GENE.in(INDEL_ALB, INDEL_SFTPB, INDEL_SLC34A2, MUTATION_EGFR))
                .fetch();

        for(Record record : result)
        {
            final String sampleId = record.getValue(SOMATICVARIANT.SAMPLEID);
            final String gene = record.getValue(SOMATICVARIANT.GENE);
            final VariantType type = VariantType.valueOf(record.getValue(SOMATICVARIANT.TYPE));

            int repeatCount = record.getValue(SOMATICVARIANT.REPEATCOUNT);
            int position = record.getValue(SOMATICVARIANT.POSITION);
            final String ref = record.getValue(SOMATICVARIANT.REF);
            final String alt = record.getValue(SOMATICVARIANT.ALT);

            if((checkIndels && isKnownIndel(gene, repeatCount, type)) || isKnownMutation(gene, type, position, ref, alt))
            {
                List<String> genes = sampleMutationMap.get(sampleId);

                if(genes == null)
                    sampleMutationMap.put(sampleId, Lists.newArrayList(gene));
                else if(!genes.contains(gene))
                    genes.add(gene);
            }
        }

        return sampleMutationMap;
    }

    private static List<String> loadSpecificMutations(final String vcfFile, boolean checkIndels)
    {
        final List<SomaticVariant> variants = SomaticDataLoader.loadSomaticVariantsFromVcf(vcfFile, Lists.newArrayList());

        final List<String> mutations = Lists.newArrayList();

        for(SomaticVariant variant : variants)
        {
            if(mutations.contains(variant.Gene))
                continue;

            if((checkIndels && isKnownIndel(variant.Gene, variant.RepeatCount, variant.Type))
            || isKnownMutation(variant.Gene, variant.Type, variant.Position, variant.Ref, variant.Alt))
            {
                mutations.add(variant.Gene);
            }
        }

        return mutations;
    }

    private static Map<String,List<DriverCatalog>> getAllDrivers(final DatabaseAccess dbAccess, final String specificSampleId)
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
                    .transcript(record.getValue(DRIVERCATALOG.TRANSCRIPTID))
                    .isCanonical(record.getValue(DRIVERCATALOG.CANONICALTRANSCRIPT) == 1)
                    .chromosome(record.getValue(DRIVERCATALOG.CHROMOSOME))
                    .chromosomeBand(record.getValue(DRIVERCATALOG.CHROMOSOMEBAND))
                    .driver(DriverType.checkConvertType(record.getValue(DRIVERCATALOG.DRIVER)))
                    .category(DriverCategory.valueOf(record.getValue(DRIVERCATALOG.CATEGORY)))
                    .likelihoodMethod(LikelihoodMethod.valueOf(record.getValue(DRIVERCATALOG.LIKELIHOODMETHOD)))
                    .driverLikelihood(record.getValue(DRIVERCATALOG.DRIVERLIKELIHOOD))
                    .missense(record.getValue(DRIVERCATALOG.MISSENSE))
                    .nonsense(record.getValue(DRIVERCATALOG.NONSENSE))
                    .splice(record.getValue(DRIVERCATALOG.SPLICE))
                    .inframe(record.getValue(DRIVERCATALOG.INFRAME))
                    .frameshift(record.getValue(DRIVERCATALOG.FRAMESHIFT))
                    .biallelic(record.getValue(DRIVERCATALOG.BIALLELIC) != 0)
                    .minCopyNumber(record.getValue(DRIVERCATALOG.MINCOPYNUMBER))
                    .maxCopyNumber(record.getValue(DRIVERCATALOG.MAXCOPYNUMBER))
                    .build();

            // ignore germline drivers
            if(DriverType.isGermline(driverCatalog.driver()))
                continue;

            final List<DriverCatalog> drivers = sampleDriverMap.get(sampleId);
            if(drivers == null)
                sampleDriverMap.put(sampleId, Lists.newArrayList(driverCatalog));
            else
                drivers.add(driverCatalog);
        }

        return sampleDriverMap;
    }

    private static void mapFeatureData(
            final String sampleId, final Map<String,List<SampleFeatureData>> sampleDrivers, final List<DriverCatalog> drivers,
            final List<LinxFusion> fusions, final List<AnnotatedVirus> virusAnnotations, final List<String> indelGenes)
    {
        final List<SampleFeatureData> featuresList = Lists.newArrayList();

        if(drivers != null)
        {
            for(final DriverCatalog driver : drivers)
            {
                if(DriverType.isGermline(driver.driver()))
                    continue;

                FeatureType featureType = driver.driver() == DriverType.AMP || driver.driver() == DriverType.PARTIAL_AMP ? AMP : DRIVER;

                SampleFeatureData feature = new SampleFeatureData(sampleId, driver.gene(), featureType, driver.driverLikelihood());

                if(driver.driver() == DriverType.AMP || driver.driver() == DriverType.PARTIAL_AMP)
                {
                    feature.ExtraInfo.put(DRIVER_TYPE, DRIVER_TYPE_AMP);
                    feature.ExtraInfo.put(AMP_CN, String.format("%.1f", driver.maxCopyNumber()));
                }
                else if(driver.driver() == DriverType.DEL)
                {
                    feature.ExtraInfo.put(DRIVER_TYPE, DRIVER_TYPE_DEL);
                }

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

        if(virusAnnotations != null)
        {
            final List<SampleFeatureData> viralInsertDataList = Lists.newArrayList();

            // convert to virus group name and ensure no duplicates per sample
            for(AnnotatedVirus annotatedVirus : virusAnnotations)
            {
                if(!annotatedVirus.reported())
                    continue;

                double likelihood = (annotatedVirus.virusDriverLikelihoodType() == HIGH || annotatedVirus.virusDriverLikelihoodType() == UNKNOWN) ? 1 : 0.5;

                // String virusName = annotatedVirus.interpretation(); // to be confirmed
                String virusName = ViralInsertionType.fromVirusName(annotatedVirus.name()).toString();

                if(viralInsertDataList.stream().anyMatch(x -> x.Name.equals(virusName)))
                    continue;

                viralInsertDataList.add(new SampleFeatureData(sampleId, virusName, FeatureType.VIRUS, likelihood));
            }

            /*
            virusAnnotations.stream()
                    .map(x -> new SampleFeatureData(sampleId, x.interpretation(), FeatureType.VIRUS, 1))
                    .filter(x -> !x.Name.equals(OTHER.toString()))
                    .filter(x -> viralInsertDataList.stream().noneMatch(y -> y.Name.equals(x.Name)))
                    .forEach(x -> viralInsertDataList.add(x));
            */

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

    public static boolean loadRefPrevalenceData(
            final String filename, final Map<String,FeaturePrevCounts> featurePrevTotals,
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

                FeaturePrevCounts prevTotals = featurePrevTotals.get(prevData.typeName());

                if(prevTotals == null)
                {
                    prevTotals = new FeaturePrevCounts();
                    featurePrevTotals.put(prevData.typeName(), prevTotals);
                }

                prevTotals.MaxPrevalence = max(prevTotals.MaxPrevalence, prevData.Prevalence);

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
            CUP_LOGGER.error("failed to read feature prevalence data: {}", e.toString());
            return false;
        }

        return true;
    }

    public static List<FeaturePrevData> loadRefFeatureOverrides(final String filename)
    {
        final List<FeaturePrevData> cancerFeaturePrevOverrides = Lists.newArrayList();

        if(filename == null || filename.isEmpty())
            return cancerFeaturePrevOverrides;

        try
        {
            final List<String> fileData = Files.readAllLines(new File(filename).toPath());

            fileData.remove(0);

            for(final String line : fileData)
            {
                final FeaturePrevData prevData = FeaturePrevData.from(line);

                if(prevData != null)
                {
                    cancerFeaturePrevOverrides.add(prevData);
                }
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read feature overrides prevalence data file: {}", e.toString());
        }

        return cancerFeaturePrevOverrides;
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

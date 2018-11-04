//package com.hartwig.hmftools.actionability;
//
//import java.io.File;
//import java.io.FileNotFoundException;
//import java.io.IOException;
//import java.nio.file.Files;
//import java.nio.file.Path;
//import java.util.List;
//import java.util.Set;
//import java.util.stream.Collectors;
//
//import com.hartwig.hmftools.common.actionability.cancertype.CancerTypeAnalyzer;
//import com.hartwig.hmftools.common.actionability.cancertype.PrimaryTumorToDOIDMapping;
//import com.hartwig.hmftools.common.actionability.cnv.CopyNumberEvidenceAnalyzer;
//import com.hartwig.hmftools.common.actionability.somaticvariant.SomaticVariantEvidenceAnalyzer;
//import com.hartwig.hmftools.common.context.ProductionRunContextFactory;
//import com.hartwig.hmftools.common.context.RunContext;
//import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
//import com.hartwig.hmftools.common.io.path.PathExtensionFinder;
//import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
//import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
//import com.hartwig.hmftools.common.hgvsCodingImpact.SomaticVariant;
//import com.hartwig.hmftools.common.hgvsCodingImpact.SomaticVariantFactory;
//
//import org.apache.commons.cli.CommandLine;
//import org.apache.commons.cli.CommandLineParser;
//import org.apache.commons.cli.DefaultParser;
//import org.apache.commons.cli.Options;
//import org.apache.commons.cli.ParseException;
//import org.apache.logging.log4j.LogManager;
//import org.apache.logging.log4j.Logger;
//import org.apache.logging.log4j.util.Strings;
//import org.jetbrains.annotations.NotNull;
//import org.jetbrains.annotations.Nullable;
//
//public class ActionabilityApplication {
//
//    private static final Logger LOGGER = LogManager.getLogger(ActionabilityApplication.class);
//    private static final String TUMOR_LOCATION_CSV = "tumor_location_csv";
//    private static final String RUN_DIRECTORY = "run_dir";
//    private static final String SOMATIC_VCF_EXTENSION_V3 = "_post_processed_v2.2.vcf.gz";
//    private static final String SOMATIC_VCF_EXTENSION_V4 = "_post_processed.vcf.gz";
//    private static final String PURPLE_DIRECTORY = "purple";
//
//    public static void main(final String... args) throws ParseException, IOException {
//        LOGGER.info("Determining actionability variants");
//        final Options options = createOptions();
//        final CommandLine cmd = createCommandLine(options, args);
//        final String runDir = cmd.getOptionValue(RUN_DIRECTORY);
//
//        final RunContext run = ProductionRunContextFactory.fromRunDirectory(runDir);
//        final List<SomaticVariant> variants = loadPassedSomaticVariants(run.tumorSample(), runDir);
//        final List<GeneCopyNumber> exomeGeneCopyNumbers = loadPurpleGeneCopyNumbers(runDir, run.tumorSample());
//        final String patientIdentifier = toPatientIdentifier(run.tumorSample());
//        LOGGER.info("Tumor sample: " + run.tumorSample());
//        LOGGER.info("patientId: " + patientIdentifier);
//
//        LOGGER.info("DOIDs of tumor location");
//        final List<PatientTumorLocation> patientTumorLocations = PatientTumorLocation.readRecords(cmd.getOptionValue(TUMOR_LOCATION_CSV));
//        final PatientTumorLocation patientTumorLocation = extractPatientTumorLocation(patientTumorLocations, run.tumorSample());
//        final String primaryTumorLocation = patientTumorLocation != null ? patientTumorLocation.primaryTumorLocation() : Strings.EMPTY;
//        final String cancerSubtype = patientTumorLocation != null ? patientTumorLocation.cancerSubtype() : Strings.EMPTY;
//        LOGGER.info("Tumor location from patient: " + primaryTumorLocation);
//        LOGGER.info("Cancer subtype from patient: " + cancerSubtype);
//
//        PrimaryTumorToDOIDMapping cancerTypeMappingReading = PrimaryTumorToDOIDMapping.createFromResource();
//        String doidsPrimaryTumorLocation = cancerTypeMappingReading.findDoids(primaryTumorLocation);
//        LOGGER.info("DOID: " + doidsPrimaryTumorLocation);
//
//        String fileCancerTumorsWithDOID = "/data/common/dbs/knowledgebases/output/knowledgebaseCancerTypes.tsv";
//
//        LOGGER.info("");
//        LOGGER.info("Start processing actionability somaticVariants");
//
//        String fileActionabilityVariants = "/data/common/dbs/knowledgebases/output/actionableVariants.tsv";
//        String fileActionabilityRanges = "/data/common/dbs/knowledgebases/output/actionableRanges.tsv";
//
//        LOGGER.info("Variants: " + variants.size());
//        if (Files.exists(new File(fileActionabilityVariants).toPath()) && Files.exists(new File(fileActionabilityRanges).toPath())
//                && Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
//            SomaticVariantEvidenceAnalyzer analyzer =
//                    SomaticVariantEvidenceAnalyzer.loadFromFileVariantsAndFileRanges(fileActionabilityVariants, fileActionabilityRanges);
//            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.createFromKnowledgeBase(fileCancerTumorsWithDOID);
//
//            Set<String> actionableGenesVariants = analyzer.actionableGenes();
//            LOGGER.info(actionableGenesVariants.size() + " actionable genes found for variants (and ranges)");
//            LOGGER.info("");
//
//            List<SomaticVariant> variantsOnActionableGenes =
//                    variants.stream().filter(hgvsCodingImpact -> actionableGenesVariants.contains(hgvsCodingImpact.gene())).collect(Collectors.toList());
//
//            LOGGER.info(
//                    "Gene" + "\t" + "Chromosome" + "\t" + "Position" + "\t" + "Ref" + "\t" + "Alt" + "\t" + "Source" + "\t" + "Drug" + "\t"
//                            + "DrugsType" + "\t" + "CancerType" + "\t" + "Level" + "\t" + "Response" + "\t" + "Actionable_variant");
//            for (SomaticVariant hgvsCodingImpact : variantsOnActionableGenes) {
//                analyzer.actionableVariants(hgvsCodingImpact, cancerTypeAnalyzer, doidsPrimaryTumorLocation);
//            }
//
//            LOGGER.info("Gene" + "\t" + "Chromosome" + "\t" + "Start" + "\t" + "End" + "\t" + " Source" + "\t" + "Drug" + "\t" + "DrugType"
//                    + "\t" + "CancerType" + "\t" + "Level" + "\t" + "Response" + "\t" + "Actionable_variant");
//            for (SomaticVariant hgvsCodingImpact : variantsOnActionableGenes) {
//                analyzer.actionableRange(hgvsCodingImpact, cancerTypeAnalyzer, doidsPrimaryTumorLocation);
//            }
//
//        } else if (!Files.exists(new File(fileActionabilityVariants).toPath())) {
//            LOGGER.warn("File does not exist: " + fileActionabilityVariants);
//        } else if (!Files.exists(new File(fileActionabilityRanges).toPath())) {
//            LOGGER.warn("File does not exist: " + fileActionabilityRanges);
//        } else if (!Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
//            LOGGER.warn("File does not exist: " + fileCancerTumorsWithDOID);
//        }
//
//        LOGGER.info("Finish processing actionability somaticVariants");
//        LOGGER.info("");
//        LOGGER.info("Start processing actionability cnvs");
//        String fileActionabilityCNVs = "/data/common/dbs/knowledgebases/output/evidenceForCopyNumber.tsv";
//
//        LOGGER.info("CNVs: " + exomeGeneCopyNumbers.size());
//        if (Files.exists(new File(fileActionabilityCNVs).toPath()) && Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
//            CopyNumberEvidenceAnalyzer analyzerCNVs = CopyNumberEvidenceAnalyzer.loadFromFileCNVs(fileActionabilityCNVs);
//            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.createFromKnowledgeBase(fileCancerTumorsWithDOID);
//
//            Set<String> actionableGenesCNVS = analyzerCNVs.actionableGenes();
//            LOGGER.info(actionableGenesCNVS.size() + " actionable genes found for cnvs");
//            LOGGER.info("");
//
//            List<GeneCopyNumber> variantsOnActionableGenes = exomeGeneCopyNumbers.stream()
//                    .filter(geneCopyNumber -> actionableGenesCNVS.contains(geneCopyNumber.gene()))
//                    .collect(Collectors.toList());
//
//            LOGGER.info(
//                    "Gene" + "\t" + "CnvType" + "\t" + "Source" + "\t" + "Drug" + "\t" + "DrugType" + "\t" + "CancerType" + "\t" + "Level"
//                            + "\t" + "Response" + "\t" + "Actionable_variant");
//            for (GeneCopyNumber geneCopyNumber : variantsOnActionableGenes) {
//                analyzerCNVs.evidenceForCopyNumber(geneCopyNumber, cancerTypeAnalyzer, doidsPrimaryTumorLocation);
//            }
//        } else if (!Files.exists(new File(fileActionabilityCNVs).toPath())) {
//            LOGGER.warn("File does not exist: " + fileActionabilityCNVs);
//        } else if (!Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
//            LOGGER.warn("File does not exist: " + fileCancerTumorsWithDOID);
//        }
//        LOGGER.info("Finished processing actionability cnvs");
//
//        //        LOGGER.info("");
//        //        LOGGER.info("Start processing fusions");
//        //        String fileActionabilityFusionPairs = "/data/common/dbs/knowledgebases/output/actionableFusionPairs.tsv";
//        //        String fileActionabilityPromiscuousFive = "/data/common/dbs/knowledgebases/output/actionablePromiscuousFive.tsv";
//        //        String fileActionabilityPromiscuousThree = "/data/common/dbs/knowledgebases/output/actionablePromiscuousThree.tsv";
//        //        if (Files.exists(new File(fileActionabilityFusionPairs).toPath()) && Files.exists(new File(fileActionabilityPromiscuousFive).toPath()) &&
//        //                Files.exists(new File(fileActionabilityPromiscuousThree).toPath()) && Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
//        //
//        //
//        //            FusionEvidenceAnalyzer analyzerFusion = FusionEvidenceAnalyzer.loadFromFileFusions(fileActionabilityFusionPairs,
//        //                    fileActionabilityPromiscuousFive, fileActionabilityPromiscuousThree);
//        //            CancerTypeAnalyzer cancerTypeAnalyzer = CancerTypeAnalyzer.createFromKnowledgeBase(fileCancerTumorsWithDOID);
//        //
//        //
//        //        } else if (!Files.exists(new File(fileActionabilityFusionPairs).toPath())) {
//        //            LOGGER.warn("File does not exist: " + fileActionabilityFusionPairs);
//        //        } else if (!Files.exists(new File(fileActionabilityPromiscuousFive).toPath())) {
//        //            LOGGER.warn("File does not exist: " + fileActionabilityPromiscuousFive);
//        //        } else if (!Files.exists(new File(fileActionabilityPromiscuousThree).toPath())) {
//        //            LOGGER.warn("File does not exist: " + fileActionabilityPromiscuousThree);
//        //        } else if (!Files.exists(new File(fileCancerTumorsWithDOID).toPath())) {
//        //            LOGGER.warn("File does not exist: " + fileCancerTumorsWithDOID);
//        //        }
//        //
//        //
//        //        LOGGER.info("Finished processing fusions");
//    }
//
//    @NotNull
//    private static List<GeneCopyNumber> loadPurpleGeneCopyNumbers(@NotNull final String runDirectory, @NotNull final String sample)
//            throws IOException {
//        final String cnvBasePath = runDirectory + File.separator + PURPLE_DIRECTORY;
//        final String fileName = GeneCopyNumberFile.generateFilename(cnvBasePath, sample);
//        return GeneCopyNumberFile.read(fileName);
//    }
//
//    @NotNull
//    private static Options createOptions() {
//        final Options options = new Options();
//        options.addOption(RUN_DIRECTORY, true, "Complete path towards a single run dir where patient reporter will run on.");
//        options.addOption(TUMOR_LOCATION_CSV, true, "Complete path towards the (curated) tumor location csv.");
//        return options;
//    }
//
//    @NotNull
//    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
//        final CommandLineParser parser = new DefaultParser();
//        return parser.parse(options, args);
//    }
//
//    @NotNull
//    private static List<SomaticVariant> loadPassedSomaticVariants(@NotNull final String sample, @NotNull final String path)
//            throws IOException {
//        // TODO (KODU): Clean up once pipeline v3 no longer exists
//        Path vcfPath;
//        try {
//            vcfPath = PathExtensionFinder.build().findPath(path, SOMATIC_VCF_EXTENSION_V3);
//        } catch (FileNotFoundException exception) {
//            vcfPath = PathExtensionFinder.build().findPath(path, SOMATIC_VCF_EXTENSION_V4);
//        }
//        return SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, vcfPath.toString());
//    }
//
//    @Nullable
//    private static PatientTumorLocation extractPatientTumorLocation(@NotNull final List<PatientTumorLocation> patientTumorLocations,
//            @NotNull final String sample) {
//        final String patientIdentifier = toPatientIdentifier(sample);
//
//        final List<PatientTumorLocation> matchingIdTumorLocations = patientTumorLocations.stream()
//                .filter(patientTumorLocation -> patientTumorLocation.patientIdentifier().equals(patientIdentifier))
//                .collect(Collectors.toList());
//
//        // KODU: We should never have more than one curated tumor location for a single patient.
//        assert matchingIdTumorLocations.size() < 2;
//
//        if (matchingIdTumorLocations.size() == 1) {
//            return matchingIdTumorLocations.get(0);
//        } else {
//            LOGGER.warn("Could not find patient " + patientIdentifier + " in clinical data!");
//            return null;
//        }
//    }
//
//    @NotNull
//    private static String toPatientIdentifier(@NotNull final String sample) {
//        if (sample.length() >= 12 && (sample.startsWith("CPCT") || sample.startsWith("DRUP"))) {
//            return sample.substring(0, 12);
//        }
//        // KODU: If we want to generate a report for non-CPCT/non-DRUP we assume patient and sample are identical.
//        return sample;
//    }
//}

package com.hartwig.hmftools.ckb;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.ckb.clinicaltrials.ClinicalTrial;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CkbImporterApplication {

    private static final Logger LOGGER = LogManager.getLogger(CkbImporterApplication.class);
    public static final String VERSION = CkbImporterApplication.class.getPackage().getImplementationVersion();

    private static final String CLINICAL_TRIALS = "clinicalTrials";
    private static final String DRUG_CLASSES = "drugClasses";
    private static final String DRUGS = "drugs";
    private static final String GENES = "genes";
    private static final String GLOBAL_THERAPY_APPROVAL_STATUSES = "globalTherapyApprovalStatuses";
    private static final String INDICATIONS = "indications";
    private static final String MOLECULAR_PROFILES = "molecularProfiles";
    private static final String REFERENCES = "references";
    private static final String THERAPIES = "therapies";
    private static final String TREATMENT_APPROACHES = "treatmentApproaches";
    private static final String VARIANTS = "variants";

    public static void main(String[] args) throws IOException {
        LOGGER.info("Running CKB importer v{}", VERSION);

        Options options = CkbImporterConfig.createOptions();

        CkbImporterConfig config = null;
        try {
            config = CkbImporterConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("CKB Importer", options);
            System.exit(1);
        }

        readJsonData(config);
        LOGGER.info("Complete!");
    }

    private static void readJsonData(@NotNull CkbImporterConfig config) throws IOException {
        List<ClinicalTrial> clinicalTrials = readingClinicalTrial(config.cbkDir() + CLINICAL_TRIALS);
        readingDrugsClasses(config.cbkDir() + DRUG_CLASSES);
        readingDrugs(config.cbkDir() + DRUGS);
        readingGenes(config.cbkDir() + GENES);
        readingGlobalTherapyApprovalStatuses(config.cbkDir() + GLOBAL_THERAPY_APPROVAL_STATUSES);
        readingIndications(config.cbkDir() + INDICATIONS);
        readingMolecularProfiles(config.cbkDir() + MOLECULAR_PROFILES);
        readingReferences(config.cbkDir() + REFERENCES);
        readingTherapies(config.cbkDir() + THERAPIES);
        readingTreatmentApproaches(config.cbkDir() + TREATMENT_APPROACHES);
        readingVariants(config.cbkDir() + VARIANTS);
    }

    @NotNull
    private static List<ClinicalTrial> readingClinicalTrial(@NotNull String clinicalTrialDir) throws IOException {
        LOGGER.info("Start reading clinical trials");

        List<ClinicalTrial> clinicalTrials = Lists.newArrayList();
        File[] filesClinicalTrials = new File(clinicalTrialDir).listFiles();
        if (filesClinicalTrials != null) {
            LOGGER.info("The total files in the clinical trial dir is {}", filesClinicalTrials.length);
            for (File clinicalTrial : filesClinicalTrials) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(clinicalTrial));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject clinicalTrialsEntryObject = parser.parse(reader).getAsJsonObject();
                    //clinicalTrials.add();
                }
            }
        }
        LOGGER.info("Finished reading clinical trials");

        return clinicalTrials;
    }

    private static void readingDrugsClasses(@NotNull String drugsClassesDir) throws IOException {
        LOGGER.info("Start reading drugs classes");

        LOGGER.info("Finished reading drugs classes");
    }

    private static void readingDrugs(@NotNull String drugsDir) throws IOException {
        LOGGER.info("Start reading drugs");

        LOGGER.info("Finished reading drugs");
    }

    private static void readingGenes(@NotNull String genesDir) throws IOException {
        LOGGER.info("Start reading genes");

        LOGGER.info("Finished reading genes");
    }

    private static void readingGlobalTherapyApprovalStatuses(@NotNull String globalTherpyApprovalStatusesDir) throws IOException {
        LOGGER.info("Start reading global therpy approval statuses");

        LOGGER.info("Finished reading global therpy approval statuses");
    }

    private static void readingIndications(@NotNull String indicationsDir) throws IOException {
        LOGGER.info("Start reading indications");

        LOGGER.info("Finished reading indications");
    }

    private static void readingMolecularProfiles(@NotNull String molecularProfilesDir) throws IOException {
        LOGGER.info("Start reading molecular profiles");

        LOGGER.info("Finished reading molecular profiles");
    }

    private static void readingReferences(@NotNull String referencesDir) throws IOException {
        LOGGER.info("Start reading references");

        LOGGER.info("Finished reading references");
    }

    private static void readingTreatmentApproaches(@NotNull String therapiesDir) throws IOException {
        LOGGER.info("Start reading therapies");

        LOGGER.info("Finished reading therapies");
    }

    private static void readingTherapies(@NotNull String treatmentApproachesDir) throws IOException {
        LOGGER.info("Start reading treatment approaches");

        LOGGER.info("Finished reading treatment approaches");
    }

    private static void readingVariants(@NotNull String variantsDir) throws IOException {
        LOGGER.info("Start reading variants");

        LOGGER.info("Finished reading variants");
    }
}

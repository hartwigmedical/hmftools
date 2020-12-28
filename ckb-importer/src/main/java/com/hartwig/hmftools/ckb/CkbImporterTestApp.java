package com.hartwig.hmftools.ckb;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.net.InetAddress;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.ckb.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.clinicaltrial.ClinicalTrialLocations;
import com.hartwig.hmftools.ckb.clinicaltrial.ImmutableClinicalTrial;
import com.hartwig.hmftools.ckb.clinicaltrial.ImmutableClinicalTrialLocations;
import com.hartwig.hmftools.ckb.clinicaltrial.ImmutableIndications;
import com.hartwig.hmftools.ckb.clinicaltrial.ImmutableTherapies;
import com.hartwig.hmftools.ckb.clinicaltrial.Indications;
import com.hartwig.hmftools.ckb.clinicaltrial.Therapies;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class CkbImporterTestApp {

    private static final Logger LOGGER = LogManager.getLogger(CkbImporterTestApp.class);

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
        Configurator.setRootLevel(Level.DEBUG);

        String hostname = InetAddress.getLocalHost().getHostName();
        LOGGER.debug("Running on '{}'", hostname);

        String ckbPath;

        if (hostname.toLowerCase().contains("datastore")) {
            ckbPath = "/data/common/dbs/ckb/";
        } else {
            ckbPath = System.getProperty("user.home") + "/hmf/projects/serve/ckb/";
        }

        readJsonData(ckbPath);
        LOGGER.info("Complete!");

    }

    private static void readJsonData(@NotNull String ckbPath) throws IOException {
        List<ClinicalTrial> clinicalTrials = readingClinicalTrial(ckbPath + CLINICAL_TRIALS);
        readingDrugsClasses(ckbPath + DRUG_CLASSES);
        readingDrugs(ckbPath + DRUGS);
        readingGenes(ckbPath + GENES);
        readingGlobalTherapyApprovalStatuses(ckbPath + GLOBAL_THERAPY_APPROVAL_STATUSES);
        readingIndications(ckbPath + INDICATIONS);
        readingMolecularProfiles(ckbPath + MOLECULAR_PROFILES);
        readingReferences(ckbPath + REFERENCES);
        readingTherapies(ckbPath + THERAPIES);
        readingTreatmentApproaches(ckbPath + TREATMENT_APPROACHES);
        readingVariants(ckbPath + VARIANTS);
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
                    clinicalTrials.add(ImmutableClinicalTrial.builder()
                            .nctId(clinicalTrialsEntryObject.getAsJsonPrimitive("nctId").getAsString())
                            .title(clinicalTrialsEntryObject.getAsJsonPrimitive("title").getAsString())
                            .phase(clinicalTrialsEntryObject.getAsJsonPrimitive("phase").getAsString())
                            .recruitment(clinicalTrialsEntryObject.getAsJsonPrimitive("recruitment").getAsString())
                            .therapies(clinicalTrialsEntryObject.has("therapies") ? retrieveClinicalTrialsTherpaies(
                                    clinicalTrialsEntryObject.getAsJsonArray("therapies")) : null)
                            .ageGroups(Lists.newArrayList())
                            .gender(Strings.EMPTY)
                            .variantRequirements(clinicalTrialsEntryObject.getAsJsonPrimitive("variantRequirements").getAsString())
                            .sponsors(Strings.EMPTY)
                            .updateDate(clinicalTrialsEntryObject.getAsJsonPrimitive("updateDate").getAsString())
                            .indications(retrieveClinicalTrialsIndications(clinicalTrialsEntryObject.getAsJsonArray("indications")))
                            .variantRequirementDetails(Lists.newArrayList())
                            .clinicalTrialLocations(retrieveClinicalTrialsLocations(clinicalTrialsEntryObject.getAsJsonArray(
                                    "clinicalTrialLocations")))
                            .build());
                }
            }
        } LOGGER.info("Finished reading clinical trials");

        return clinicalTrials;
    }

    private static List<Therapies> retrieveClinicalTrialsTherpaies(@NotNull JsonArray jsonArray) {
        List<Therapies> therapies = Lists.newArrayList();
        for (JsonElement therapy : jsonArray) {
            therapies.add(ImmutableTherapies.builder()
                    .id(therapy.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .therapyName(therapy.getAsJsonObject().getAsJsonPrimitive("therapyName").getAsString())
                    .synonyms(Strings.EMPTY)
                    .build());
        }
        return therapies;
    }

    private static List<Indications> retrieveClinicalTrialsIndications(@NotNull JsonArray jsonArray) {
        List<Indications> indications = Lists.newArrayList();
        for (JsonElement indication : jsonArray) {
            indications.add(ImmutableIndications.builder()
                    .id(indication.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(indication.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .source(indication.getAsJsonObject().getAsJsonPrimitive("source").getAsString())
                    .build());
        }
        return indications;
    }

    private static List<ClinicalTrialLocations> retrieveClinicalTrialsLocations(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialLocations> locations = Lists.newArrayList();

        for (JsonElement location : jsonArray) {
            locations.add(ImmutableClinicalTrialLocations.builder()
                    .nctId(location.getAsJsonObject().getAsJsonPrimitive("nctId").getAsString())
                    .facility(Strings.EMPTY)
                    .city(location.getAsJsonObject().getAsJsonPrimitive("city").getAsString())
                    .country(location.getAsJsonObject().getAsJsonPrimitive("country").getAsString())
                    .status(Strings.EMPTY)
                    .state(Strings.EMPTY)
                    .zip(Strings.EMPTY)
                    .clinicalTrialContacts(Lists.newArrayList())

                    .build());
        }
        return locations;
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

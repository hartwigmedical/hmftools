package com.hartwig.hmftools.ckb;

import java.io.IOException;
import java.net.InetAddress;
import java.util.List;

import com.hartwig.hmftools.ckb.clinicaltrials.ClinicalTrial;
import com.hartwig.hmftools.ckb.clinicaltrials.ClinicalTrialFactory;
import com.hartwig.hmftools.ckb.drugclasses.DrugClass;
import com.hartwig.hmftools.ckb.drugclasses.DrugClassFactory;
import com.hartwig.hmftools.ckb.drugs.Drugs;
import com.hartwig.hmftools.ckb.drugs.DrugsFactory;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

public class CkbImporterTestApp {

    private static final Logger LOGGER = LogManager.getLogger(CkbImporterTestApp.class);

    private static final String CLINICAL_TRIALS = "clinicalTrials";
    private static final String DRUG_CLASSES = "drugclasses";
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
        List<ClinicalTrial> clinicalTrials = ClinicalTrialFactory.readingClinicalTrial(ckbPath + CLINICAL_TRIALS);
        List<DrugClass> drugClasses = DrugClassFactory.readingDrugClasses(ckbPath + DRUG_CLASSES);
        List<Drugs> drugs = DrugsFactory.readingDrugs(ckbPath + DRUGS);
//        readingGenes(ckbPath + GENES);
//        readingGlobalTherapyApprovalStatuses(ckbPath + GLOBAL_THERAPY_APPROVAL_STATUSES);
//        readingIndications(ckbPath + INDICATIONS);
//        readingMolecularProfiles(ckbPath + MOLECULAR_PROFILES);
//        readingReferences(ckbPath + REFERENCES);
//        readingTherapies(ckbPath + THERAPIES);
//        readingTreatmentApproaches(ckbPath + TREATMENT_APPROACHES);
//        readingVariants(ckbPath + VARIANTS);
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

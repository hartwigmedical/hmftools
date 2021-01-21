package com.hartwig.hmftools.ckb;

import java.io.IOException;
import java.net.InetAddress;
import java.util.List;

import com.hartwig.hmftools.ckb.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.clinicaltrial.ClinicalTrialFactory;
import com.hartwig.hmftools.ckb.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.drugclass.DrugClassFactory;
import com.hartwig.hmftools.ckb.drug.Drug;
import com.hartwig.hmftools.ckb.drug.DrugFactory;
import com.hartwig.hmftools.ckb.gene.Gene;
import com.hartwig.hmftools.ckb.gene.GeneFactory;
import com.hartwig.hmftools.ckb.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.globaltherapyapprovalstatus.GlobalTherapyApprovalStatusFactory;
import com.hartwig.hmftools.ckb.indication.Indication;
import com.hartwig.hmftools.ckb.indication.IndicationFactory;
import com.hartwig.hmftools.ckb.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.molecularprofile.MolecularprofileFactory;

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
//        List<ClinicalTrial> clinicalTrials = ClinicalTrialFactory.readingClinicalTrial(ckbPath + CLINICAL_TRIALS);
//        List<DrugClass> drugClasses = DrugClassFactory.readingDrugClasses(ckbPath + DRUG_CLASSES);
//        List<Drug> drugs = DrugFactory.readingDrugs(ckbPath + DRUGS);
//        List<Gene> genes = GeneFactory.readingGenes(ckbPath + GENES);
//        List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatuses =
//                GlobalTherapyApprovalStatusFactory.readingGlobalTherapyApprovalStatus(ckbPath + GLOBAL_THERAPY_APPROVAL_STATUSES);
//        List<Indication> indications = IndicationFactory.readingIndication(ckbPath + INDICATIONS);
         List<MolecularProfile> molecularProfiles = MolecularprofileFactory.readingMolecularprofile(ckbPath + MOLECULAR_PROFILES);
        //        readingReferences(ckbPath + REFERENCES);
        //        readingTherapies(ckbPath + THERAPIES);
        //        readingTreatmentApproaches(ckbPath + TREATMENT_APPROACHES);
        //        readingVariants(ckbPath + VARIANTS);
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

package com.hartwig.hmftools.ckb;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.ckb.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.clinicaltrial.ClinicalTrialFactory;
import com.hartwig.hmftools.ckb.drug.Drug;
import com.hartwig.hmftools.ckb.drug.DrugFactory;
import com.hartwig.hmftools.ckb.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.drugclass.DrugClassFactory;
import com.hartwig.hmftools.ckb.gene.Gene;
import com.hartwig.hmftools.ckb.gene.GeneFactory;
import com.hartwig.hmftools.ckb.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.globaltherapyapprovalstatus.GlobalTherapyApprovalStatusFactory;
import com.hartwig.hmftools.ckb.indication.Indication;
import com.hartwig.hmftools.ckb.indication.IndicationFactory;
import com.hartwig.hmftools.ckb.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.molecularprofile.MolecularprofileFactory;
import com.hartwig.hmftools.ckb.reference.Reference;
import com.hartwig.hmftools.ckb.reference.ReferenceFactory;
import com.hartwig.hmftools.ckb.therapy.Therapy;
import com.hartwig.hmftools.ckb.therapy.TherapyFactory;
import com.hartwig.hmftools.ckb.treatmentApproach.TreatmentApproach;
import com.hartwig.hmftools.ckb.treatmentApproach.TreatmentApproachFactory;
import com.hartwig.hmftools.ckb.variant.Variant;
import com.hartwig.hmftools.ckb.variant.VariantFactory;

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

        CkbEntry ckbEntry = readJsonData(config);
        writingDataToDatabase(ckbEntry);

        LOGGER.info("Complete!");
    }

    @NotNull
    private static CkbEntry readJsonData(@NotNull CkbImporterConfig config) throws IOException {
        LOGGER.info("Start with reading all files");

        String ckbPath = config.cbkDir();

        List<ClinicalTrial> clinicalTrials = ClinicalTrialFactory.readingClinicalTrial(ckbPath + File.separator + CLINICAL_TRIALS);
        List<Drug> drugs = DrugFactory.readingDrugs(ckbPath + File.separator + DRUGS);
        List<DrugClass> drugClasses = DrugClassFactory.readingDrugClasses(ckbPath + File.separator + DRUG_CLASSES);
        List<Gene> genes = GeneFactory.readingGenes(ckbPath + File.separator + GENES);
        List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatuses =
                GlobalTherapyApprovalStatusFactory.readingGlobalTherapyApprovalStatus(
                        ckbPath + File.separator + GLOBAL_THERAPY_APPROVAL_STATUSES);
        List<Indication> indications = IndicationFactory.readingIndication(ckbPath + File.separator + INDICATIONS);
        List<MolecularProfile> molecularProfiles =
                MolecularprofileFactory.readingMolecularprofile(ckbPath + File.separator + MOLECULAR_PROFILES);
        List<Reference> references = ReferenceFactory.readingReference(ckbPath + File.separator + REFERENCES);
        List<Therapy> therapies = TherapyFactory.readingTherapy(ckbPath + File.separator + THERAPIES);
        List<TreatmentApproach> treatmentApproaches =
                TreatmentApproachFactory.readingTreatmentApproch(ckbPath + File.separator + TREATMENT_APPROACHES);
        List<Variant> variants = VariantFactory.readingVariant(ckbPath + File.separator + VARIANTS);

        LOGGER.info("All files are readed");

        return ImmutableCkbEntry.builder()
                .clinicalTrial(clinicalTrials)
                .drug(drugs)
                .drugClass(drugClasses)
                .gene(genes)
                .globalTherapyApprovalStatus(globalTherapyApprovalStatuses)
                .indication(indications)
                .molecularProfile(molecularProfiles)
                .reference(references)
                .therapy(therapies)
                .treatmentApproach(treatmentApproaches)
                .variant(variants)
                .build();
    }

    private static void writingDataToDatabase(@NotNull CkbEntry ckbEntry) {
        LOGGER.info("Start with writing data to DB");

        LOGGER.info("ClinicalTrial {}", ckbEntry.clinicalTrial().get(0));
        LOGGER.info("Drug {}", ckbEntry.drug().get(0));
        LOGGER.info("DrugClass {}", ckbEntry.drugClass().get(0));
        LOGGER.info("Gene {}", ckbEntry.gene().get(0));
        LOGGER.info("GlobalTherapyApprovalStatus {}", ckbEntry.globalTherapyApprovalStatus().get(0));
        LOGGER.info("Indication {}", ckbEntry.indication().get(0));
        LOGGER.info("MolecularProfile {}", ckbEntry.molecularProfile().get(0));
        LOGGER.info("Reference {}", ckbEntry.reference().get(0));
        LOGGER.info("Therapy {}", ckbEntry.therapy().get(0));
        LOGGER.info("TreatmentApproachInfo {}", ckbEntry.treatmentApproach().get(0));
        LOGGER.info("Variant {}", ckbEntry.variant().get(0));

        LOGGER.info("All data is written to the DB");
    }

}

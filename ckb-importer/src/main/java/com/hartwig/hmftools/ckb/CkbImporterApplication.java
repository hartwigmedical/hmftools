package com.hartwig.hmftools.ckb;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.reader.clinicaltrial.ClinicalTrialFactory;
import com.hartwig.hmftools.ckb.dao.CkbDAO;
import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.reader.drug.DrugFactory;
import com.hartwig.hmftools.ckb.datamodel.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.reader.drugclass.DrugClassFactory;
import com.hartwig.hmftools.ckb.datamodel.gene.Gene;
import com.hartwig.hmftools.ckb.reader.gene.GeneFactory;
import com.hartwig.hmftools.ckb.datamodel.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.reader.globaltherapyapprovalstatus.GlobalTherapyApprovalStatusFactory;
import com.hartwig.hmftools.ckb.datamodel.indication.Indication;
import com.hartwig.hmftools.ckb.reader.indication.IndicationFactory;
import com.hartwig.hmftools.ckb.datamodel.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.reader.molecularprofile.MolecularprofileFactory;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.reader.reference.ReferenceFactory;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;
import com.hartwig.hmftools.ckb.reader.therapy.TherapyFactory;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproach.TreatmentApproach;
import com.hartwig.hmftools.ckb.reader.treatmentapproch.TreatmentApproachFactory;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.ckb.reader.variant.VariantFactory;

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

    public static void main(String[] args) throws IOException, SQLException {
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

        if (config.skipDatabaseWriting()) {
            LOGGER.info("Skipping DB writing.");
        } else {
            CkbDAO ckbDAO = connect(config);
            LOGGER.info("Deleting all from CKB db");
            ckbDAO.deleteAll();
            LOGGER.info("Starting insertion of all CKB entries");
            ckbDAO.writeCkb(ckbEntry);

        }

        LOGGER.info("Complete!");
    }

    private static CkbDAO connect(@NotNull CkbImporterConfig config) throws SQLException {
        return CkbDAO.connectToCkbDAO(config.dbUser(), config.dbPass(), "jdbc:" + config.dbUrl());
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

}

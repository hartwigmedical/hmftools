package com.hartwig.hmftools.ckb.json;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialReader;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrial;
import com.hartwig.hmftools.ckb.json.drug.DrugReader;
import com.hartwig.hmftools.ckb.json.drug.JsonDrug;
import com.hartwig.hmftools.ckb.json.drugclass.DrugClassReader;
import com.hartwig.hmftools.ckb.json.drugclass.JsonDrugClass;
import com.hartwig.hmftools.ckb.json.gene.GeneReader;
import com.hartwig.hmftools.ckb.json.gene.JsonGene;
import com.hartwig.hmftools.ckb.json.globaltherapyapprovalstatus.GlobalTherapyApprovalStatusReader;
import com.hartwig.hmftools.ckb.json.globaltherapyapprovalstatus.JsonGlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.json.indication.IndicationReader;
import com.hartwig.hmftools.ckb.json.indication.JsonIndication;
import com.hartwig.hmftools.ckb.json.molecularprofile.JsonMolecularProfile;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfileReader;
import com.hartwig.hmftools.ckb.json.reference.JsonReference;
import com.hartwig.hmftools.ckb.json.reference.ReferenceReader;
import com.hartwig.hmftools.ckb.json.therapy.JsonTherapy;
import com.hartwig.hmftools.ckb.json.therapy.TherapyReader;
import com.hartwig.hmftools.ckb.json.treatmentapproach.JsonTreatmentApproach;
import com.hartwig.hmftools.ckb.json.treatmentapproach.TreatmentApproachReader;
import com.hartwig.hmftools.ckb.json.variant.JsonVariant;
import com.hartwig.hmftools.ckb.json.variant.VariantReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CkbJsonReader {

    private static final Logger LOGGER = LogManager.getLogger(CkbJsonReader.class);

    private static final String MOLECULAR_PROFILES_DIR = "molecularProfiles";
    private static final String VARIANTS_DIR = "variants";
    private static final String GENES_DIR = "genes";
    private static final String INDICATIONS_DIR = "indications";
    private static final String TREATMENT_APPROACHES_DIR = "treatmentApproaches";
    private static final String THERAPIES_DIR = "therapies";
    private static final String DRUGS_DIR = "drugs";
    private static final String DRUG_CLASSES_DIR = "drugClasses";
    private static final String CLINICAL_TRIALS_DIR = "clinicalTrials";
    private static final String GLOBAL_THERAPY_APPROVAL_STATUSES_DIR = "globalTherapyApprovalStatuses";
    private static final String REFERENCES_DIR = "references";

    private CkbJsonReader() {
    }

    @NotNull
    public static CkbJsonDatabase read(@NotNull String ckbDir) throws IOException {
        return read(ckbDir, null);
    }

    @NotNull
    public static CkbJsonDatabase read(@NotNull String ckbDir, @Nullable Integer maxFilesToReadPerType) throws IOException {
        List<JsonMolecularProfile> molecularProfiles =
                new MolecularProfileReader(maxFilesToReadPerType).read(ckbDir + File.separator + MOLECULAR_PROFILES_DIR);
        List<JsonVariant> variants = new VariantReader(maxFilesToReadPerType).read(ckbDir + File.separator + VARIANTS_DIR);
        List<JsonGene> genes = new GeneReader(maxFilesToReadPerType).read(ckbDir + File.separator + GENES_DIR);
        List<JsonIndication> indications = new IndicationReader(maxFilesToReadPerType).read(ckbDir + File.separator + INDICATIONS_DIR);
        List<JsonTreatmentApproach> treatmentApproaches =
                new TreatmentApproachReader(maxFilesToReadPerType).read(ckbDir + File.separator + TREATMENT_APPROACHES_DIR);
        List<JsonTherapy> therapies = new TherapyReader(maxFilesToReadPerType).read(ckbDir + File.separator + THERAPIES_DIR);
        List<JsonDrug> drugs = new DrugReader(maxFilesToReadPerType).read(ckbDir + File.separator + DRUGS_DIR);
        List<JsonDrugClass> drugClasses = new DrugClassReader(maxFilesToReadPerType).read(ckbDir + File.separator + DRUG_CLASSES_DIR);
        List<JsonClinicalTrial> clinicalTrials =
                new ClinicalTrialReader(maxFilesToReadPerType).read(ckbDir + File.separator + CLINICAL_TRIALS_DIR);
        List<JsonGlobalTherapyApprovalStatus> globalTherapyApprovalStatuses =
                new GlobalTherapyApprovalStatusReader(maxFilesToReadPerType).read(
                        ckbDir + File.separator + GLOBAL_THERAPY_APPROVAL_STATUSES_DIR);
        List<JsonReference> references = new ReferenceReader(maxFilesToReadPerType).read(ckbDir + File.separator + REFERENCES_DIR);

        return ImmutableCkbJsonDatabase.builder()
                .molecularProfiles(molecularProfiles)
                .variants(variants)
                .genes(genes)
                .indications(indications)
                .treatmentApproaches(treatmentApproaches)
                .therapies(therapies)
                .drugs(drugs)
                .drugClasses(drugClasses)
                .clinicalTrials(clinicalTrials)
                .globalTherapyApprovalStatuses(globalTherapyApprovalStatuses)
                .references(references)
                .build();
    }
}

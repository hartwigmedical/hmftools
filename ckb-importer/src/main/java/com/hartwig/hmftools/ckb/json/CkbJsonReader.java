package com.hartwig.hmftools.ckb.json;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialReader;
import com.hartwig.hmftools.ckb.json.drug.Drug;
import com.hartwig.hmftools.ckb.json.drug.DrugReader;
import com.hartwig.hmftools.ckb.json.drugclass.DrugClass;
import com.hartwig.hmftools.ckb.json.drugclass.DrugClassReader;
import com.hartwig.hmftools.ckb.json.gene.Gene;
import com.hartwig.hmftools.ckb.json.gene.GeneReader;
import com.hartwig.hmftools.ckb.json.globaltherapyapprovalstatus.GlobalTherapyApprovalStatus;
import com.hartwig.hmftools.ckb.json.globaltherapyapprovalstatus.GlobalTherapyApprovalStatusReader;
import com.hartwig.hmftools.ckb.json.indication.Indication;
import com.hartwig.hmftools.ckb.json.indication.IndicationReader;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfileReader;
import com.hartwig.hmftools.ckb.json.reference.Reference;
import com.hartwig.hmftools.ckb.json.reference.ReferenceReader;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;
import com.hartwig.hmftools.ckb.json.therapy.TherapyReader;
import com.hartwig.hmftools.ckb.json.treatmentapproch.TreatmentApproach;
import com.hartwig.hmftools.ckb.json.treatmentapproch.TreatmentApproachReader;
import com.hartwig.hmftools.ckb.json.variant.Variant;
import com.hartwig.hmftools.ckb.json.variant.VariantReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CkbJsonReader {

    private static final Logger LOGGER = LogManager.getLogger(CkbJsonReader.class);

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

    private CkbJsonReader() {
    }

    @NotNull
    public static CkbJsonDatabase read(@NotNull String ckbDir) throws IOException {
        return read(ckbDir, null);
    }

    @NotNull
    public static CkbJsonDatabase read(@NotNull String ckbDir, @Nullable Integer maxFilesToReadPerType) throws IOException {
        LOGGER.info("Reading all CKB json files from '{}'", ckbDir);

        List<ClinicalTrial> clinicalTrials = new ClinicalTrialReader(maxFilesToReadPerType).read(ckbDir + File.separator + CLINICAL_TRIALS);
        List<Drug> drugs = new DrugReader(maxFilesToReadPerType).read(ckbDir + File.separator + DRUGS);
        List<DrugClass> drugClasses = new DrugClassReader(maxFilesToReadPerType).read(ckbDir + File.separator + DRUG_CLASSES);
        List<Gene> genes = new GeneReader(maxFilesToReadPerType).read(ckbDir + File.separator + GENES);
        List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatuses = new GlobalTherapyApprovalStatusReader(maxFilesToReadPerType).read(
                ckbDir + File.separator + GLOBAL_THERAPY_APPROVAL_STATUSES);
        List<Indication> indications = new IndicationReader(maxFilesToReadPerType).read(ckbDir + File.separator + INDICATIONS);
        List<MolecularProfile> molecularProfiles =
                new MolecularProfileReader(maxFilesToReadPerType).read(ckbDir + File.separator + MOLECULAR_PROFILES);
        List<Reference> references = new ReferenceReader(maxFilesToReadPerType).read(ckbDir + File.separator + REFERENCES);
        List<Therapy> therapies = new TherapyReader(maxFilesToReadPerType).read(ckbDir + File.separator + THERAPIES);
        List<TreatmentApproach> treatmentApproaches =
                new TreatmentApproachReader(maxFilesToReadPerType).read(ckbDir + File.separator + TREATMENT_APPROACHES);
        List<Variant> variants = new VariantReader(maxFilesToReadPerType).read(ckbDir + File.separator + VARIANTS);

        LOGGER.info(" Json file reading completed");

        return ImmutableCkbJsonDatabase.builder()
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

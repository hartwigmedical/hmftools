package com.hartwig.hmftools.ckb;

import static com.hartwig.hmftools.ckb.database.Tables.DRUG;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGCLINICALTRIAL;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGCLINICALTRIALTHERAPY;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGDESCRIPTION;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGDESCRIPTIONREFERENCE;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGDRUGCLASS;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGEVIDENCE;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGEVIDENCEINDICATION;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGEVIDENCEMOLECULARPROFILE;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGEVIDENCEREFERENCE;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGEVIDENCETHERAPY;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGEVIDENCETREATMENTAPPROCH;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGGLOBALAPPROAVALSTATUS;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGGLOBALAPPROAVALSTATUSINDICATION;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGGLOBALAPPROAVALSTATUSMOLECULARPROFILE;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGGLOBALAPPROAVALSTATUSTHERAPY;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGSYNONYM;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGTERM;
import static com.hartwig.hmftools.ckb.database.Tables.DRUGTHERAPY;

import java.util.List;

import com.hartwig.hmftools.ckb.datamodel.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.datamodel.common.DescriptionInfo;
import com.hartwig.hmftools.ckb.datamodel.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.datamodel.common.EvidenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.datamodel.common.IndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.datamodel.drug.Drug;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;

public class DrugDAO {

    private DrugDAO() {
    }

    public static void writeDrug(@NotNull DSLContext context, @NotNull Drug drug) {
        int drugId =
                context.insertInto(DRUG, DRUG.IDOFDRUG, DRUG.DRUGNAME, DRUG.TRADENAME, DRUG.CASREGISTRYNUM, DRUG.NCTID, DRUG.CREATEDATE)
                        .values(drug.id(), drug.drugName(), drug.tradeName(), drug.casRegistryNum(), drug.nctId(), drug.createDate())
                        .returning(DRUG.ID)
                        .fetchOne()
                        .getValue(DRUG.ID);
        writeDrugTerm(drugId, drug.term(), context);
        writeDrugSynonyms(drugId, drug.synonym(), context);
        writeDrugDescription(drugId, drug.description(), context);
        writeDrugDrugClass(drugId, drug.drugClass(), context);
        writeDrugClinicalTrial(drugId, drug.clinicalTrial(), context);
        writeDrugEvidence(drugId, drug.evidence(), context);
        writeDrugTherapy(drugId, drug.therapy(), context);
        writeDrugGlobalApprovalStatus(drugId, drug.globalApprovalStatus(), context);

    }

    private static void writeDrugTerm(int drugId, @NotNull List<String> terms, @NotNull DSLContext context) {
        for (String term : terms) {
            context.insertInto(DRUGTERM, DRUGTERM.DRUGID, DRUGTERM.TERM).values(drugId, term).execute();
        }
    }

    private static void writeDrugSynonyms(int drugId, @NotNull List<String> synonyms, @NotNull DSLContext context) {
        for (String synonym : synonyms) {
            context.insertInto(DRUGSYNONYM, DRUGSYNONYM.DRUGID, DRUGSYNONYM.SYNONYM).values(drugId, synonym).execute();
        }
    }

    private static void writeDrugDescription(int drugId, @NotNull List<DescriptionInfo> descriptions, @NotNull DSLContext context) {
        for (DescriptionInfo descriptionInfo : descriptions) {
            int desciptionId = context.insertInto(DRUGDESCRIPTION, DRUGDESCRIPTION.DRUGID, DRUGDESCRIPTION.DESCRIPTION)
                    .values(drugId, descriptionInfo.description())
                    .returning(DRUGDESCRIPTION.ID)
                    .fetchOne()
                    .getValue(DRUGDESCRIPTION.ID);
            writeDescriptionReference(desciptionId, descriptionInfo.reference(), context);
        }
    }

    private static void writeDescriptionReference(int descriptionId, @NotNull List<ReferenceInfo> references, @NotNull DSLContext context) {
        for (ReferenceInfo reference : references) {
            context.insertInto(DRUGDESCRIPTIONREFERENCE,
                    DRUGDESCRIPTIONREFERENCE.DRUGDESCRIPTIONID,
                    DRUGDESCRIPTIONREFERENCE.DRUGDESCRIPTIONREFERENCEID,
                    DRUGDESCRIPTIONREFERENCE.PUBMEDID,
                    DRUGDESCRIPTIONREFERENCE.TITLE,
                    DRUGDESCRIPTIONREFERENCE.URL,
                    DRUGDESCRIPTIONREFERENCE.AUTHORS,
                    DRUGDESCRIPTIONREFERENCE.JOURNAL,
                    DRUGDESCRIPTIONREFERENCE.VOLUME,
                    DRUGDESCRIPTIONREFERENCE.ISSUE,
                    DRUGDESCRIPTIONREFERENCE.DATE,
                    DRUGDESCRIPTIONREFERENCE.ABSTRACTTEXT,
                    DRUGDESCRIPTIONREFERENCE.YEAR)
                    .values(descriptionId,
                            reference.id(),
                            reference.pubMedId(),
                            reference.title(),
                            reference.url(),
                            reference.authors(),
                            reference.journal(),
                            reference.volume(),
                            reference.issue(),
                            reference.date(),
                            reference.abstractText(),
                            reference.year())
                    .execute();
        }
    }

    private static void writeDrugDrugClass(int drugId, @NotNull List<DrugClassInfo> drugClasses, @NotNull DSLContext context) {
        for (DrugClassInfo drugClass : drugClasses) {
            context.insertInto(DRUGDRUGCLASS, DRUGDRUGCLASS.DRUGID, DRUGDRUGCLASS.DRUGDRUGCLASSID, DRUGDRUGCLASS.DRUGCLASS)
                    .values(drugId, drugClass.id(), drugClass.drugClass())
                    .execute();
        }
    }

    private static void writeDrugClinicalTrial(int drugId, @NotNull List<ClinicalTrialInfo> drugClinicalTrials,
            @NotNull DSLContext context) {
        for (ClinicalTrialInfo drugClinicalTrial : drugClinicalTrials) {
            int drugClinicalTrialId = context.insertInto(DRUGCLINICALTRIAL,
                    DRUGCLINICALTRIAL.DRUGID,
                    DRUGCLINICALTRIAL.NCTID,
                    DRUGCLINICALTRIAL.TITLE,
                    DRUGCLINICALTRIAL.PHASE,
                    DRUGCLINICALTRIAL.RECRUITMENT)
                    .values(drugId,
                            drugClinicalTrial.nctId(),
                            drugClinicalTrial.title(),
                            drugClinicalTrial.phase(),
                            drugClinicalTrial.recruitment())
                    .returning(DRUGCLINICALTRIAL.ID)
                    .fetchOne()
                    .getValue(DRUGCLINICALTRIAL.ID);
            writeDrugClinicalTrialTherapy(drugClinicalTrialId, drugClinicalTrial.therapy(), context);
        }
    }

    private static void writeDrugClinicalTrialTherapy(int drugClinicalTrialId, @NotNull List<TherapyInfo> clinicalTrailtherapies,
            @NotNull DSLContext context) {
        for (TherapyInfo drugClinicalTrialTherapy : clinicalTrailtherapies) {
            context.insertInto(DRUGCLINICALTRIALTHERAPY,
                    DRUGCLINICALTRIALTHERAPY.DRUGCLINICALTRIALID,
                    DRUGCLINICALTRIALTHERAPY.DRUGCLINICALTRIALTHERAPYID,
                    DRUGCLINICALTRIALTHERAPY.THERAPYNAME,
                    DRUGCLINICALTRIALTHERAPY.SYNONYM)
                    .values(drugClinicalTrialId,
                            drugClinicalTrialTherapy.id(),
                            drugClinicalTrialTherapy.therapyName(),
                            drugClinicalTrialTherapy.synonyms())
                    .execute();
        }
    }

    private static void writeDrugEvidence(int drugId, @NotNull List<EvidenceInfo> drugEvidences, @NotNull DSLContext context) {
        for (EvidenceInfo drugEvidence : drugEvidences) {
            int drugEvidenceId = context.insertInto(DRUGEVIDENCE,
                    DRUGEVIDENCE.DRUGID,
                    DRUGEVIDENCE.DRUGEVIDENCEID,
                    DRUGEVIDENCE.APPROVALSTATUS,
                    DRUGEVIDENCE.EVIDENCETYPE,
                    DRUGEVIDENCE.EFFICACYEVIDENCE,
                    DRUGEVIDENCE.RESPONSETYPE,
                    DRUGEVIDENCE.AMPCAPASCOEVIDENCELEVEL,
                    DRUGEVIDENCE.AMPCAPASCOINFERREDTIER)
                    .values(drugId,
                            drugEvidence.id(),
                            drugEvidence.approvalStatus(),
                            drugEvidence.evidenceType(),
                            drugEvidence.efficacyEvidence(),
                            drugEvidence.responseType(),
                            drugEvidence.ampCapAscoEvidenceLevel(),
                            drugEvidence.ampCapAscoInferredTier())
                    .returning(DRUGEVIDENCE.ID)
                    .fetchOne()
                    .getValue(DRUGEVIDENCE.ID);
            writeDrugEvidenceMolecularProfile(drugEvidenceId, drugEvidence.molecularProfile(), context);
            writeDrugEvidenceTherapy(drugEvidenceId, drugEvidence.therapy(), context);
            writeDrugEvidenceIndication(drugEvidenceId, drugEvidence.indication(), context);
            writeDrugEvidenceReference(drugEvidenceId, drugEvidence.reference(), context);
            writeDrugEvidenceTreatmentApproach(drugEvidenceId, drugEvidence.treatmentApproach(), context);
        }
    }

    private static void writeDrugEvidenceMolecularProfile(int drugEvidenceId, @NotNull MolecularProfileInfo drugEvidenceMolecularProfile,
            @NotNull DSLContext context) {
        context.insertInto(DRUGEVIDENCEMOLECULARPROFILE,
                DRUGEVIDENCEMOLECULARPROFILE.DRUGEVIDENCEID,
                DRUGEVIDENCEMOLECULARPROFILE.DRUGEVIDENCEMOLECULARPROFILEID,
                DRUGEVIDENCEMOLECULARPROFILE.PROFILENAME)
                .values(drugEvidenceId, drugEvidenceMolecularProfile.id(), drugEvidenceMolecularProfile.profileName())
                .execute();
    }

    private static void writeDrugEvidenceTherapy(int drugEvidenceId, @NotNull TherapyInfo drugEvidenceTherapy,
            @NotNull DSLContext context) {
        context.insertInto(DRUGEVIDENCETHERAPY,
                DRUGEVIDENCETHERAPY.DRUGEVIDENCEID,
                DRUGEVIDENCETHERAPY.DRUGEVIDENCETHERAPYID,
                DRUGEVIDENCETHERAPY.THERAPYNAME,
                DRUGEVIDENCETHERAPY.SYNONYM)
                .values(drugEvidenceId, drugEvidenceTherapy.id(), drugEvidenceTherapy.therapyName(), drugEvidenceTherapy.synonyms())
                .execute();
    }

    private static void writeDrugEvidenceIndication(int drugEvidenceId, @Nullable IndicationInfo drugEvidenceIndication,
            @NotNull DSLContext context) {
        if (drugEvidenceIndication != null) {
            context.insertInto(DRUGEVIDENCEINDICATION,
                    DRUGEVIDENCEINDICATION.DRUGEVIDENCEID,
                    DRUGEVIDENCEINDICATION.DRUGEVIDENCEINDICATIONID,
                    DRUGEVIDENCEINDICATION.NAME,
                    DRUGEVIDENCEINDICATION.SOURCE)
                    .values(drugEvidenceId, drugEvidenceIndication.id(), drugEvidenceIndication.name(), drugEvidenceIndication.source())
                    .execute();
        }
    }

    private static void writeDrugEvidenceReference(int drugEvidenceId, @NotNull List<ReferenceInfo> drugEvidenceReferences,
            @NotNull DSLContext context) {
        for (ReferenceInfo drugEvidenceReference : drugEvidenceReferences) {
            context.insertInto(DRUGEVIDENCEREFERENCE,
                    DRUGEVIDENCEREFERENCE.DRUGEVIDENCEID,
                    DRUGEVIDENCEREFERENCE.DRUGEVIDENCEREFERENCEID,
                    DRUGEVIDENCEREFERENCE.PUBMEDID,
                    DRUGEVIDENCEREFERENCE.TITLE,
                    DRUGEVIDENCEREFERENCE.URL)
                    .values(drugEvidenceId,
                            drugEvidenceReference.id(),
                            drugEvidenceReference.pubMedId(),
                            drugEvidenceReference.title(),
                            drugEvidenceReference.url())
                    .execute();
        }
    }

    private static void writeDrugEvidenceTreatmentApproach(int drugEvidenceId,
            @Nullable List<TreatmentApproachInfo> drugEvidenceTreatmentApproaches, @NotNull DSLContext context) {
        if (drugEvidenceTreatmentApproaches != null) {
            for (TreatmentApproachInfo drugEvidenceTreatmentApproach : drugEvidenceTreatmentApproaches) {
                context.insertInto(DRUGEVIDENCETREATMENTAPPROCH,
                        DRUGEVIDENCETREATMENTAPPROCH.DRUGEVIDENCEID,
                        DRUGEVIDENCETREATMENTAPPROCH.DRUGEVIDENCETREATMENTAPPROCHID,
                        DRUGEVIDENCETREATMENTAPPROCH.NAME,
                        DRUGEVIDENCETREATMENTAPPROCH.PROFILENAME)
                        .values(drugEvidenceId,
                                drugEvidenceTreatmentApproach.id(),
                                drugEvidenceTreatmentApproach.name(),
                                drugEvidenceTreatmentApproach.profileName())
                        .execute();
            }
        }
    }

    private static void writeDrugTherapy(int drugId, @NotNull List<TherapyInfo> drugTherapies, @NotNull DSLContext context) {
        for (TherapyInfo drugTherapy : drugTherapies) {
            context.insertInto(DRUGTHERAPY, DRUGTHERAPY.DRUGID, DRUGTHERAPY.DRUGTHERAPYID, DRUGTHERAPY.THERAPYNAME, DRUGTHERAPY.SYNONYM)
                    .values(drugId, drugTherapy.id(), drugTherapy.therapyName(), drugTherapy.synonyms())
                    .execute();
        }
    }

    private static void writeDrugGlobalApprovalStatus(int drugId, @Nullable List<GlobalApprovalStatusInfo> drugGlobalApproavalStatusses,
            @NotNull DSLContext context) {
        if (drugGlobalApproavalStatusses != null) {
            for (GlobalApprovalStatusInfo drugGlobalApprovalStatus : drugGlobalApproavalStatusses) {
                int drugGlobalApprovalStatusId = context.insertInto(DRUGGLOBALAPPROAVALSTATUS,
                        DRUGGLOBALAPPROAVALSTATUS.DRUGID,
                        DRUGGLOBALAPPROAVALSTATUS.DRUGGLOBALAPPROAVALSTATUSID,
                        DRUGGLOBALAPPROAVALSTATUS.APPROVALAUTHORITY,
                        DRUGGLOBALAPPROAVALSTATUS.APPROVALSTATUS)
                        .values(drugId,
                                drugGlobalApprovalStatus.id(),
                                drugGlobalApprovalStatus.approvalAuthority(),
                                drugGlobalApprovalStatus.approvalStatus())
                        .returning(DRUGGLOBALAPPROAVALSTATUS.ID)
                        .fetchOne()
                        .getValue(DRUGGLOBALAPPROAVALSTATUS.ID);
                writeDrugGlobalApproavalStatusTherapy(drugGlobalApprovalStatusId, drugGlobalApprovalStatus.therapy(), context);
                writeDrugGlobalApproavalStatusIndication(drugGlobalApprovalStatusId, drugGlobalApprovalStatus.indication(), context);
                writeDrugGlobalApprovalStatusMolecularProfile(drugGlobalApprovalStatusId,
                        drugGlobalApprovalStatus.molecularProfile(),
                        context);
            }
        }
    }

    private static void writeDrugGlobalApproavalStatusTherapy(int drugGlopbalApproavalStatusId,
            @NotNull TherapyInfo drugGlobalApproavalStatusTherapy, @NotNull DSLContext context) {
        context.insertInto(DRUGGLOBALAPPROAVALSTATUSTHERAPY,
                DRUGGLOBALAPPROAVALSTATUSTHERAPY.DRUGGLOBALAPPROAVALSTATUSID,
                DRUGGLOBALAPPROAVALSTATUSTHERAPY.DRUGGLOBALAPPROAVALSTATUSTHERAPYID,
                DRUGGLOBALAPPROAVALSTATUSTHERAPY.THERAPYNAME,
                DRUGGLOBALAPPROAVALSTATUSTHERAPY.SYNONYM)
                .values(drugGlopbalApproavalStatusId,
                        drugGlobalApproavalStatusTherapy.id(),
                        drugGlobalApproavalStatusTherapy.therapyName(),
                        drugGlobalApproavalStatusTherapy.synonyms())
                .execute();
    }

    private static void writeDrugGlobalApproavalStatusIndication(int drugGlobalApprovalStatusId,
            @NotNull IndicationInfo drugGlobalApprovalStatusIndication, @NotNull DSLContext context) {
        context.insertInto(DRUGGLOBALAPPROAVALSTATUSINDICATION,
                DRUGGLOBALAPPROAVALSTATUSINDICATION.DRUGGLOBALAPPROAVALSTATUSID,
                DRUGGLOBALAPPROAVALSTATUSINDICATION.DRUGGLOBALAPPROAVALSTATUSINDICATIONID,
                DRUGGLOBALAPPROAVALSTATUSINDICATION.NAME,
                DRUGGLOBALAPPROAVALSTATUSINDICATION.SOURCE)
                .values(drugGlobalApprovalStatusId,
                        drugGlobalApprovalStatusIndication.id(),
                        drugGlobalApprovalStatusIndication.name(),
                        drugGlobalApprovalStatusIndication.source())
                .execute();
    }

    private static void writeDrugGlobalApprovalStatusMolecularProfile(int drugGlobalApprovalStatusId,
            @NotNull MolecularProfileInfo drugGlobalApprovalStatusMolecularProfile, @NotNull DSLContext context) {
        context.insertInto(DRUGGLOBALAPPROAVALSTATUSMOLECULARPROFILE,
                DRUGGLOBALAPPROAVALSTATUSMOLECULARPROFILE.DRUGGLOBALAPPROAVALSTATUSID,
                DRUGGLOBALAPPROAVALSTATUSMOLECULARPROFILE.DRUGGLOBALAPPROAVALSTATUSMOLECULARPROFILEID,
                DRUGGLOBALAPPROAVALSTATUSMOLECULARPROFILE.PROFILENAME)
                .values(drugGlobalApprovalStatusId,
                        drugGlobalApprovalStatusMolecularProfile.id(),
                        drugGlobalApprovalStatusMolecularProfile.profileName())
                .execute();
    }

    public static void clearDrug(@NotNull DSLContext context) {
        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
        context.truncate(DRUGGLOBALAPPROAVALSTATUSMOLECULARPROFILE).execute();
        context.truncate(DRUGGLOBALAPPROAVALSTATUSINDICATION).execute();
        context.truncate(DRUGGLOBALAPPROAVALSTATUSTHERAPY).execute();
        context.truncate(DRUGGLOBALAPPROAVALSTATUS).execute();
        context.truncate(DRUGTHERAPY).execute();
        context.truncate(DRUGEVIDENCETREATMENTAPPROCH).execute();
        context.truncate(DRUGEVIDENCEREFERENCE).execute();
        context.truncate(DRUGEVIDENCEINDICATION).execute();
        context.truncate(DRUGEVIDENCETHERAPY).execute();
        context.truncate(DRUGEVIDENCEMOLECULARPROFILE).execute();
        context.truncate(DRUGEVIDENCE).execute();
        context.truncate(DRUGCLINICALTRIALTHERAPY).execute();
        context.truncate(DRUGCLINICALTRIAL).execute();
        context.truncate(DRUGDRUGCLASS).execute();
        context.truncate(DRUGDESCRIPTIONREFERENCE).execute();
        context.truncate(DRUGDESCRIPTION).execute();
        context.truncate(DRUGSYNONYM).execute();
        context.truncate(DRUGTERM).execute();
        context.truncate(DRUG).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }
}

package com.hartwig.hmftools.ckb.dao;

import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialContact;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialLocation;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.common.IndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;

import static com.hartwig.hmftools.ckb.database.Tables.CLINICALTRIAL;
import static com.hartwig.hmftools.ckb.database.Tables.CLINICALTRIALAGEGROUP;
import static com.hartwig.hmftools.ckb.database.Tables.CLINICALTRIALCONTACT;
import static com.hartwig.hmftools.ckb.database.Tables.CLINICALTRIALINDICATION;
import static com.hartwig.hmftools.ckb.database.Tables.CLINICALTRIALLOCATION;
import static com.hartwig.hmftools.ckb.database.Tables.CLINICALTRIALTHERAPY;
import static com.hartwig.hmftools.ckb.database.Tables.CLINICALTRIALVARIANTREQUIREMENTDETAIL;
import static com.hartwig.hmftools.ckb.database.Tables.CLINICALTRIALVARIANTREQUIREMENTDETAILMOLECULARPROFILE;

import java.util.List;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class ClinicalTrialDAO {

    private ClinicalTrialDAO() {
    }

//    static void writeClinicalTrial(@NotNull DSLContext context, @NotNull ClinicalTrial clinicalTrial) {
//        int clinicalTrialId = context.insertInto(CLINICALTRIAL,
//                CLINICALTRIAL.NCTID,
//                CLINICALTRIAL.TITLE,
//                CLINICALTRIAL.PHASE,
//                CLINICALTRIAL.RECRUITMENT,
//                CLINICALTRIAL.GENDER,
//                CLINICALTRIAL.VARIANTREQUIREMENT,
//                CLINICALTRIAL.SPONSORS,
//                CLINICALTRIAL.UPDATEDATE)
//                .values(clinicalTrial.nctId(),
//                        clinicalTrial.title(),
//                        clinicalTrial.phase(),
//                        clinicalTrial.recruitment(),
//                        clinicalTrial.gender(),
//                        clinicalTrial.variantRequirement(),
//                        clinicalTrial.sponsors(),
//                        clinicalTrial.updateDate())
//                .returning(CLINICALTRIAL.ID)
//                .fetchOne()
//                .getValue(CLINICALTRIAL.ID);
//        writeAgegroups(clinicalTrialId, clinicalTrial.ageGroup(), context);
//        writeClinicalTrialTherapy(clinicalTrialId, clinicalTrial.therapy(), context);
//        writeClinicalTrialIndication(clinicalTrialId, clinicalTrial.indication(), context);
//        writeClinicalTrialVariantRequirementDetail(clinicalTrialId, clinicalTrial.variantRequirementDetail(), context);
//        writeClinicalTrialLocation(clinicalTrialId, clinicalTrial.clinicalTrialLocation(), context);
//    }
//
//    private static void writeAgegroups(int clinicalTrialId, @NotNull List<String> ageGroups, @NotNull DSLContext context) {
//        for (String ageGroup : ageGroups) {
//            context.insertInto(CLINICALTRIALAGEGROUP, CLINICALTRIALAGEGROUP.CLINICALTRIALID, CLINICALTRIALAGEGROUP.AGEGROUP)
//                    .values(clinicalTrialId, ageGroup)
//                    .execute();
//        }
//    }
//
//    private static void writeClinicalTrialTherapy(int clinicalTrialId, @NotNull List<TherapyInfo> therapies, @NotNull DSLContext context) {
//        for (TherapyInfo therapy : therapies) {
//            context.insertInto(CLINICALTRIALTHERAPY,
//                    CLINICALTRIALTHERAPY.CLINICALTRIALID,
//                    CLINICALTRIALTHERAPY.THERAPYID,
//                    CLINICALTRIALTHERAPY.THERAPYNAME,
//                    CLINICALTRIALTHERAPY.SYNONYMS)
//                    .values(clinicalTrialId, therapy.id(), therapy.therapyName(), therapy.synonyms())
//                    .execute();
//        }
//    }
//
//    private static void writeClinicalTrialIndication(int clinicalTrialId, @NotNull List<IndicationInfo> indications,
//            @NotNull DSLContext context) {
//        for (IndicationInfo indication : indications) {
//            context.insertInto(CLINICALTRIALINDICATION,
//                    CLINICALTRIALINDICATION.CLINICALTRIALID,
//                    CLINICALTRIALINDICATION.INDICATIONID,
//                    CLINICALTRIALINDICATION.NAME,
//                    CLINICALTRIALINDICATION.SOURCE)
//                    .values(clinicalTrialId, indication.id(), indication.name(), indication.source())
//                    .execute();
//        }
//    }
//
//    private static void writeClinicalTrialVariantRequirementDetail(int clinicalTrialId,
//            @NotNull List<ClinicalTrialVariantRequirementDetail> variantRequirementDetails, @NotNull DSLContext context) {
//        for (ClinicalTrialVariantRequirementDetail variantRequirementDetail : variantRequirementDetails) {
//            int variantRequirementDetailId = context.insertInto(CLINICALTRIALVARIANTREQUIREMENTDETAIL,
//                    CLINICALTRIALVARIANTREQUIREMENTDETAIL.CLINICALTRIALID,
//                    CLINICALTRIALVARIANTREQUIREMENTDETAIL.REQUIREMENTTYPE)
//                    .values(clinicalTrialId, variantRequirementDetail.requirementType())
//                    .returning(CLINICALTRIALVARIANTREQUIREMENTDETAIL.ID)
//                    .fetchOne()
//                    .getValue(CLINICALTRIALVARIANTREQUIREMENTDETAIL.ID);
//            writeClinicalTrialVariantRequirementDetailMolecularProfile(variantRequirementDetailId,
//                    variantRequirementDetail.molecularProfile(),
//                    context);
//        }
//    }
//
//    private static void writeClinicalTrialVariantRequirementDetailMolecularProfile(int variantRequirementDetailId,
//            @NotNull MolecularProfileInfo variantRequirementDetailMolecularProfile, @NotNull DSLContext context) {
//        context.insertInto(CLINICALTRIALVARIANTREQUIREMENTDETAILMOLECULARPROFILE,
//                CLINICALTRIALVARIANTREQUIREMENTDETAILMOLECULARPROFILE.CLINICALTRIALVARIANTREQUIREMENTDETAILID,
//                CLINICALTRIALVARIANTREQUIREMENTDETAILMOLECULARPROFILE.IDMOLECULARPROFILE,
//                CLINICALTRIALVARIANTREQUIREMENTDETAILMOLECULARPROFILE.PROFILENAME)
//                .values(variantRequirementDetailId,
//                        variantRequirementDetailMolecularProfile.id(),
//                        variantRequirementDetailMolecularProfile.profileName())
//                .execute();
//    }
//
//    private static void writeClinicalTrialLocation(int clinicalTrialId, @NotNull List<ClinicalTrialLocation> clinicalTrialLocations,
//            @NotNull DSLContext context) {
//        for (ClinicalTrialLocation clinicalTrialLocation : clinicalTrialLocations) {
//            int clinicalTrialLocationId = context.insertInto(CLINICALTRIALLOCATION,
//                    CLINICALTRIALLOCATION.CLINICALTRIALID,
//                    CLINICALTRIALLOCATION.NCTID,
//                    CLINICALTRIALLOCATION.FACILITY,
//                    CLINICALTRIALLOCATION.CITY,
//                    CLINICALTRIALLOCATION.COUNTRY,
//                    CLINICALTRIALLOCATION.STATUS,
//                    CLINICALTRIALLOCATION.STATE,
//                    CLINICALTRIALLOCATION.ZIP)
//                    .values(clinicalTrialId,
//                            clinicalTrialLocation.nctId(),
//                            clinicalTrialLocation.facility(),
//                            clinicalTrialLocation.city(),
//                            clinicalTrialLocation.country(),
//                            clinicalTrialLocation.status(),
//                            clinicalTrialLocation.state(),
//                            clinicalTrialLocation.zip())
//                    .returning(CLINICALTRIALLOCATION.ID)
//                    .fetchOne()
//                    .getValue(CLINICALTRIALLOCATION.ID);
//            writeClinicalTrialContact(clinicalTrialLocationId, clinicalTrialLocation.clinicalTrialContact(), context);
//        }
//    }
//
//    private static void writeClinicalTrialContact(int clinicalTrialLocationId, @NotNull List<ClinicalTrialContact> clinicalTrialContacts,
//            @NotNull DSLContext context) {
//        for (ClinicalTrialContact clinicalTrialContact : clinicalTrialContacts) {
//            context.insertInto(CLINICALTRIALCONTACT,
//                    CLINICALTRIALCONTACT.CLINICALTRIALLOCATIONID,
//                    CLINICALTRIALCONTACT.NAME,
//                    CLINICALTRIALCONTACT.EMAIL,
//                    CLINICALTRIALCONTACT.PHONE,
//                    CLINICALTRIALCONTACT.PHONEEXT,
//                    CLINICALTRIALCONTACT.ROLE)
//                    .values(clinicalTrialLocationId,
//                            clinicalTrialContact.name(),
//                            clinicalTrialContact.email(),
//                            clinicalTrialContact.phone(),
//                            clinicalTrialContact.phoneExt(),
//                            clinicalTrialContact.role())
//                    .execute();
//        }
//    }
//
//    static void clearClinicalTrial(@NotNull DSLContext context) {
//        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
//        context.truncate(CLINICALTRIALCONTACT).execute();
//        context.truncate(CLINICALTRIALLOCATION).execute();
//        context.truncate(CLINICALTRIALVARIANTREQUIREMENTDETAILMOLECULARPROFILE).execute();
//        context.truncate(CLINICALTRIALVARIANTREQUIREMENTDETAIL).execute();
//        context.truncate(CLINICALTRIALINDICATION).execute();
//        context.truncate(CLINICALTRIALTHERAPY).execute();
//        context.truncate(CLINICALTRIALAGEGROUP).execute();
//        context.truncate(CLINICALTRIAL).execute();
//        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
//    }
}

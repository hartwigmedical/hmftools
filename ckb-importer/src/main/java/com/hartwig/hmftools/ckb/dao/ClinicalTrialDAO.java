package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.tables.Agegroup.AGEGROUP;
import static com.hartwig.hmftools.ckb.database.tables.Clinicaltrial.CLINICALTRIAL;
import static com.hartwig.hmftools.ckb.database.tables.Contact.CONTACT;
import static com.hartwig.hmftools.ckb.database.tables.Indicationclinicaltrial.INDICATIONCLINICALTRIAL;
import static com.hartwig.hmftools.ckb.database.tables.Location.LOCATION;
import static com.hartwig.hmftools.ckb.database.tables.Therapyclinicaltrial.THERAPYCLINICALTRIAL;
import static com.hartwig.hmftools.ckb.database.tables.Variantrequirementdetail.VARIANTREQUIREMENTDETAIL;

import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.Contact;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.Location;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.VariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.indication.Indication;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class ClinicalTrialDAO {

    @NotNull
    private final DSLContext context;
    @NotNull
    private final TherapyDAO therapyDAO;
    @NotNull
    private final IndicationDAO indicationDAO;

    public ClinicalTrialDAO(@NotNull final DSLContext context, @NotNull final TherapyDAO therapyDAO,
            @NotNull final IndicationDAO indicationDAO) {
        this.context = context;
        this.therapyDAO = therapyDAO;
        this.indicationDAO = indicationDAO;
    }

    public void deleteAll() {
        // Note that deletions should go from branch to root
        context.deleteFrom(CONTACT).execute();
        context.deleteFrom(LOCATION).execute();

        context.deleteFrom(VARIANTREQUIREMENTDETAIL).execute();
        context.deleteFrom(AGEGROUP).execute();

        context.deleteFrom(INDICATIONCLINICALTRIAL).execute();
        context.deleteFrom(THERAPYCLINICALTRIAL).execute();

        context.deleteFrom(CLINICALTRIAL).execute();
    }

    public void write(@NotNull ClinicalTrial clinicalTrial, int ckbEntryId) {
        int id = context.insertInto(CLINICALTRIAL,
                CLINICALTRIAL.CKBENTRYID,
                CLINICALTRIAL.UPDATEDATE,
                CLINICALTRIAL.NCTID,
                CLINICALTRIAL.TITLE,
                CLINICALTRIAL.PHASE,
                CLINICALTRIAL.RECRUITMENT,
                CLINICALTRIAL.GENDER,
                CLINICALTRIAL.SPONSORS,
                CLINICALTRIAL.VARIANTREQUIREMENT)
                .values(ckbEntryId,
                        clinicalTrial.updateDate(),
                        clinicalTrial.nctId(),
                        clinicalTrial.title(),
                        clinicalTrial.phase(),
                        clinicalTrial.recruitment(),
                        clinicalTrial.gender(),
                        clinicalTrial.sponsors(),
                        clinicalTrial.variantRequirement())
                .returning(CLINICALTRIAL.ID)
                .fetchOne()
                .getValue(CLINICALTRIAL.ID);

        for (Therapy therapy : clinicalTrial.therapies()) {
            int therapyId = therapyDAO.write(therapy);
            context.insertInto(THERAPYCLINICALTRIAL, THERAPYCLINICALTRIAL.CLINICALTRIALID, THERAPYCLINICALTRIAL.THERAPYID)
                    .values(id, therapyId)
                    .execute();
        }

        for (Indication indication : clinicalTrial.indications()) {
            int indicationId = indicationDAO.write(indication);
            context.insertInto(INDICATIONCLINICALTRIAL, INDICATIONCLINICALTRIAL.CLINICALTRIALID, INDICATIONCLINICALTRIAL.INDICATIONID)
                    .values(id, indicationId)
                    .execute();
        }

        for (String ageGroup : clinicalTrial.ageGroups()) {
            context.insertInto(AGEGROUP, AGEGROUP.CLINICALTRIALID, AGEGROUP.AGEGROUP_).values(id, ageGroup).execute();
        }

        for (VariantRequirementDetail variantRequirementDetail : clinicalTrial.variantRequirementDetails()) {
            context.insertInto(VARIANTREQUIREMENTDETAIL,
                    VARIANTREQUIREMENTDETAIL.CLINICALTRIALID,
                    VARIANTREQUIREMENTDETAIL.CKBPROFILEID,
                    VARIANTREQUIREMENTDETAIL.REQUIREMENTTYPE)
                    .values(id, variantRequirementDetail.profileId(), variantRequirementDetail.requirementType())
                    .execute();
        }

        for (Location location : clinicalTrial.locations()) {
            writeLocation(location, id);
        }

    }

    private void writeLocation(@NotNull Location location, int clinicalTrialId) {
        int id = context.insertInto(LOCATION,
                LOCATION.CLINICALTRIALID,
                LOCATION.NCTID,
                LOCATION.STATUS,
                LOCATION.FACILITY,
                LOCATION.CITY,
                LOCATION.STATE,
                LOCATION.ZIP,
                LOCATION.COUNTRY)
                .values(clinicalTrialId,
                        location.nctId(),
                        location.status(),
                        location.facility(),
                        location.city(),
                        location.state(),
                        location.zip(),
                        location.country())
                .returning(LOCATION.ID)
                .fetchOne()
                .getValue(LOCATION.ID);

        for (Contact contact : location.contacts()) {
            context.insertInto(CONTACT, CONTACT.LOCATIONID, CONTACT.NAME, CONTACT.EMAIL, CONTACT.PHONE, CONTACT.PHONEEXT, CONTACT.ROLE)
                    .values(id, contact.name(), contact.email(), contact.phone(), contact.phoneExt(), contact.role())
                    .execute();
        }
    }
}

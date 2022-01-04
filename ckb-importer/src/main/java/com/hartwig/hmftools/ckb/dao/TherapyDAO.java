package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.Tables.DRUGREFERENCE;
import static com.hartwig.hmftools.ckb.database.Tables.THERAPYREFERENCE;
import static com.hartwig.hmftools.ckb.database.tables.Drug.DRUG;
import static com.hartwig.hmftools.ckb.database.tables.Drugclass.DRUGCLASS;
import static com.hartwig.hmftools.ckb.database.tables.Drugsynonym.DRUGSYNONYM;
import static com.hartwig.hmftools.ckb.database.tables.Drugterm.DRUGTERM;
import static com.hartwig.hmftools.ckb.database.tables.Globalapprovalstatus.GLOBALAPPROVALSTATUS;
import static com.hartwig.hmftools.ckb.database.tables.Therapy.THERAPY;
import static com.hartwig.hmftools.ckb.database.tables.Therapysynonym.THERAPYSYNONYM;

import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.drug.DrugClass;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.therapy.GlobalApprovalStatus;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class TherapyDAO {

    @NotNull
    private final DSLContext context;

    public TherapyDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAll() {
        // Note that deletions should go from branch to root
        context.deleteFrom(DRUGREFERENCE).execute();
        context.deleteFrom(DRUGTERM).execute();
        context.deleteFrom(DRUGSYNONYM).execute();
        context.deleteFrom(DRUGCLASS).execute();
        context.deleteFrom(DRUG).execute();

        context.deleteFrom(THERAPYSYNONYM).execute();
        context.deleteFrom(THERAPYREFERENCE).execute();

        context.deleteFrom(GLOBALAPPROVALSTATUS).execute();

        context.deleteFrom(THERAPY).execute();
    }

    public int write(@NotNull Therapy therapy) {
        int id = context.insertInto(THERAPY,
                THERAPY.CKBTHERAPYID,
                THERAPY.CREATEDATE,
                THERAPY.UPDATEDATE,
                THERAPY.THERAPYNAME,
                THERAPY.DESCRIPTION)
                .values(therapy.id(),
                        therapy.createDate(),
                        therapy.updateDate(),
                        therapy.therapyName(),
                        therapy.description())
                .returning(THERAPY.ID)
                .fetchOne()
                .getValue(THERAPY.ID);

        for (Drug drug : therapy.drugs()) {
            writeDrug(drug, id);
        }

        for (String synonym : therapy.synonyms()) {
            context.insertInto(THERAPYSYNONYM, THERAPYSYNONYM.THERAPYID, THERAPYSYNONYM.SYNONYM).values(id, synonym).execute();
        }

        for (Reference reference : therapy.references()) {
            writeTherapyReference(reference, id);
        }

        for (GlobalApprovalStatus globalApprovalStatus : therapy.globalApprovalStatuses()) {
            context.insertInto(GLOBALAPPROVALSTATUS,
                    GLOBALAPPROVALSTATUS.THERAPYID,
                    GLOBALAPPROVALSTATUS.CKBGLOBALAPPROVALSTATUSID,
                    GLOBALAPPROVALSTATUS.CKBPROFILEID,
                    GLOBALAPPROVALSTATUS.CKBINDICATIONID,
                    GLOBALAPPROVALSTATUS.APPROVALSTATUS,
                    GLOBALAPPROVALSTATUS.APPROVALAUTHORITY)
                    .values(id,
                            globalApprovalStatus.id(),
                            globalApprovalStatus.profileId(),
                            globalApprovalStatus.indicationId(),
                            globalApprovalStatus.approvalStatus(),
                            globalApprovalStatus.approvalAuthority())
                    .execute();
        }
        return id;
    }

    private void writeDrug(@NotNull Drug drug, int therapyId) {
        int id = context.insertInto(DRUG,
                DRUG.THERAPYID,
                DRUG.CKBDRUGID,
                DRUG.CREATEDATE,
                DRUG.DRUGNAME,
                DRUG.TRADENAME,
                DRUG.CASREGISTRYNUM,
                DRUG.NCITID,
                DRUG.DESCRIPTION)
                .values(therapyId,
                        drug.id(),
                        drug.createDate(),
                        drug.drugName(),
                        drug.tradeName(),
                        drug.casRegistryNum(),
                        drug.ncitId(),
                        drug.description())
                .returning(DRUG.ID)
                .fetchOne()
                .getValue(DRUG.ID);

        for (DrugClass drugClass : drug.drugClasses()) {
            context.insertInto(DRUGCLASS, DRUGCLASS.DRUGID, DRUGCLASS.CKBDRUGCLASSID, DRUGCLASS.CREATEDATE, DRUGCLASS.DRUGCLASS_)
                    .values(id, drugClass.id(), drugClass.createDate(), drugClass.drugClass())
                    .execute();
        }

        for (String term : drug.terms()) {
            context.insertInto(DRUGTERM, DRUGTERM.DRUGID, DRUGTERM.TERM).values(id, term).execute();
        }

        for (String synonym : drug.synonyms()) {
            context.insertInto(DRUGSYNONYM, DRUGSYNONYM.DRUGID, DRUGSYNONYM.SYNONYM).values(id, synonym).execute();
        }

        for (Reference reference : drug.references()) {
            writeDrugReference(reference, id);
        }
    }

    private void writeDrugReference(@NotNull Reference reference, int drugId) {
        context.insertInto(DRUGREFERENCE,
                DRUGREFERENCE.DRUGID,
                DRUGREFERENCE.CKBREFERENCEID,
                DRUGREFERENCE.PUBMEDID,
                DRUGREFERENCE.TITLE,
                DRUGREFERENCE.ABSTRACTTEXT,
                DRUGREFERENCE.URL,
                DRUGREFERENCE.JOURNAL,
                DRUGREFERENCE.AUTHORS,
                DRUGREFERENCE.VOLUME,
                DRUGREFERENCE.ISSUE,
                DRUGREFERENCE.DATE,
                DRUGREFERENCE.YEAR)
                .values(drugId,
                        reference.id(),
                        reference.pubMedId(),
                        reference.title(),
                        reference.abstractText(),
                        reference.url(),
                        reference.journal(),
                        reference.authors(),
                        reference.volume(),
                        reference.issue(),
                        reference.date(),
                        reference.year())
                .execute();
    }

    private void writeTherapyReference(@NotNull Reference reference, int therapyId) {
        context.insertInto(THERAPYREFERENCE,
                THERAPYREFERENCE.THERAPYID,
                THERAPYREFERENCE.CKBREFERENCEID,
                THERAPYREFERENCE.PUBMEDID,
                THERAPYREFERENCE.TITLE,
                THERAPYREFERENCE.ABSTRACTTEXT,
                THERAPYREFERENCE.URL,
                THERAPYREFERENCE.JOURNAL,
                THERAPYREFERENCE.AUTHORS,
                THERAPYREFERENCE.VOLUME,
                THERAPYREFERENCE.ISSUE,
                THERAPYREFERENCE.DATE,
                THERAPYREFERENCE.YEAR)
                .values(therapyId,
                        reference.id(),
                        reference.pubMedId(),
                        reference.title(),
                        reference.abstractText(),
                        reference.url(),
                        reference.journal(),
                        reference.authors(),
                        reference.volume(),
                        reference.issue(),
                        reference.date(),
                        reference.year())
                .execute();
    }
}

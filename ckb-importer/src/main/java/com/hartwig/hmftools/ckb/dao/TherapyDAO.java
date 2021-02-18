package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.tables.Drug.DRUG;
import static com.hartwig.hmftools.ckb.database.tables.Drugclass.DRUGCLASS;
import static com.hartwig.hmftools.ckb.database.tables.Drugdescription.DRUGDESCRIPTION;
import static com.hartwig.hmftools.ckb.database.tables.Drugdescriptionreference.DRUGDESCRIPTIONREFERENCE;
import static com.hartwig.hmftools.ckb.database.tables.Drugsynonym.DRUGSYNONYM;
import static com.hartwig.hmftools.ckb.database.tables.Drugterm.DRUGTERM;
import static com.hartwig.hmftools.ckb.database.tables.Globalapprovalstatus.GLOBALAPPROVALSTATUS;
import static com.hartwig.hmftools.ckb.database.tables.Therapy.THERAPY;
import static com.hartwig.hmftools.ckb.database.tables.Therapydescription.THERAPYDESCRIPTION;
import static com.hartwig.hmftools.ckb.database.tables.Therapydescriptionreference.THERAPYDESCRIPTIONREFERENCE;
import static com.hartwig.hmftools.ckb.database.tables.Therapysynonym.THERAPYSYNONYM;

import com.hartwig.hmftools.ckb.datamodel.drug.Drug;
import com.hartwig.hmftools.ckb.datamodel.drug.DrugClass;
import com.hartwig.hmftools.ckb.datamodel.drug.DrugDescription;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.therapy.GlobalApprovalStatus;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;
import com.hartwig.hmftools.ckb.datamodel.therapy.TherapyDescription;

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
        context.deleteFrom(DRUGDESCRIPTIONREFERENCE).execute();
        context.deleteFrom(DRUGDESCRIPTION).execute();
        context.deleteFrom(DRUGTERM).execute();
        context.deleteFrom(DRUGSYNONYM).execute();
        context.deleteFrom(DRUGCLASS).execute();
        context.deleteFrom(DRUG).execute();

        context.deleteFrom(THERAPYSYNONYM).execute();
        context.deleteFrom(THERAPYDESCRIPTIONREFERENCE).execute();
        context.deleteFrom(THERAPYDESCRIPTION).execute();

        context.deleteFrom(GLOBALAPPROVALSTATUS).execute();

        context.deleteFrom(THERAPY).execute();
    }

    public int write(@NotNull Therapy therapy) {
        int id = context.insertInto(THERAPY, THERAPY.CKBTHERAPYID, THERAPY.CREATEDATE, THERAPY.UPDATEDATE, THERAPY.THERAPYNAME)
                .values(therapy.id(), Util.sqlDate(therapy.createDate()), Util.sqlDate(therapy.updateDate()), therapy.therapyName())
                .returning(THERAPY.ID)
                .fetchOne()
                .getValue(THERAPY.ID);

        for (Drug drug : therapy.drugs()) {
            writeDrug(drug, id);
        }

        for (String synonym : therapy.synonyms()) {
            context.insertInto(THERAPYSYNONYM, THERAPYSYNONYM.ID, THERAPYSYNONYM.SYNONYM).values(id, synonym).execute();
        }

        for (TherapyDescription therapyDescription : therapy.descriptions()) {
            writeTherapyDescription(therapyDescription, id);
        }

        for (GlobalApprovalStatus globalApprovalStatus : therapy.globalApprovalStatuses()) {
            context.insertInto(GLOBALAPPROVALSTATUS,
                    GLOBALAPPROVALSTATUS.THERAPYID,
                    GLOBALAPPROVALSTATUS.CKBGLOBALAPPROVALSTATUSID,
                    GLOBALAPPROVALSTATUS.CKBPROFILEID,
                    GLOBALAPPROVALSTATUS.CKBTHERAPYID,
                    GLOBALAPPROVALSTATUS.CKBINDICATIONID,
                    GLOBALAPPROVALSTATUS.APPROVALSTATUS,
                    GLOBALAPPROVALSTATUS.APPROVALAUTHORITY)
                    .values(id,
                            globalApprovalStatus.id(),
                            globalApprovalStatus.profileId(),
                            globalApprovalStatus.therapyId(),
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
                DRUG.NCITID)
                .values(therapyId,
                        drug.id(),
                        Util.sqlDate(drug.createDate()),
                        drug.drugName(),
                        drug.tradeName(),
                        drug.casRegistryNum(),
                        drug.ncitId())
                .returning(DRUG.ID)
                .fetchOne()
                .getValue(DRUG.ID);

        for (DrugClass drugClass : drug.drugClasses()) {
            context.insertInto(DRUGCLASS, DRUGCLASS.DRUGID, DRUGCLASS.CKBDRUGCLASSID, DRUGCLASS.CREATEDATE, DRUGCLASS.DRUGCLASS_)
                    .values(id, drugClass.id(), Util.sqlDate(drugClass.createDate()), drugClass.drugClass())
                    .execute();
        }

        for (String term : drug.terms()) {
            context.insertInto(DRUGTERM, DRUGTERM.DRUGID, DRUGTERM.TERM).values(id, term).execute();
        }

        for (String synonym : drug.synonyms()) {
            context.insertInto(DRUGSYNONYM, DRUGSYNONYM.DRUGID, DRUGSYNONYM.SYNONYM).values(id, synonym).execute();
        }

        for (DrugDescription drugDescription : drug.descriptions()) {
            writeDrugDescription(drugDescription, id);
        }
    }

    private void writeDrugDescription(@NotNull DrugDescription drugDescription, int drugId) {
        int id = context.insertInto(DRUGDESCRIPTION, DRUGDESCRIPTION.DRUGID, DRUGDESCRIPTION.DESCRIPTION)
                .values(drugId, drugDescription.description())
                .returning(DRUGDESCRIPTION.ID)
                .fetchOne()
                .getValue(DRUGDESCRIPTION.ID);

        for (Reference reference : drugDescription.references()) {
            context.insertInto(DRUGDESCRIPTIONREFERENCE,
                    DRUGDESCRIPTIONREFERENCE.DRUGDESCRIPTIONID,
                    DRUGDESCRIPTIONREFERENCE.CKBREFERENCEID,
                    DRUGDESCRIPTIONREFERENCE.PUBMEDID,
                    DRUGDESCRIPTIONREFERENCE.TITLE,
                    DRUGDESCRIPTIONREFERENCE.ABSTRACTTEXT,
                    DRUGDESCRIPTIONREFERENCE.URL,
                    DRUGDESCRIPTIONREFERENCE.JOURNAL,
                    DRUGDESCRIPTIONREFERENCE.AUTHORS,
                    DRUGDESCRIPTIONREFERENCE.VOLUME,
                    DRUGDESCRIPTIONREFERENCE.ISSUE,
                    DRUGDESCRIPTIONREFERENCE.DATE,
                    DRUGDESCRIPTIONREFERENCE.YEAR)
                    .values(id,
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

    private void writeTherapyDescription(@NotNull TherapyDescription therapyDescription, int therapyId) {
        int id = context.insertInto(THERAPYDESCRIPTION, THERAPYDESCRIPTION.THERAPYID, THERAPYDESCRIPTION.DESCRIPTION)
                .values(therapyId, therapyDescription.description())
                .returning(THERAPYDESCRIPTION.ID)
                .fetchOne()
                .getValue(THERAPYDESCRIPTION.ID);

        for (Reference reference : therapyDescription.references()) {
            context.insertInto(THERAPYDESCRIPTIONREFERENCE,
                    THERAPYDESCRIPTIONREFERENCE.THERAPYDESCRIPTIONID,
                    THERAPYDESCRIPTIONREFERENCE.CKBREFERENCEID,
                    THERAPYDESCRIPTIONREFERENCE.PUBMEDID,
                    THERAPYDESCRIPTIONREFERENCE.TITLE,
                    THERAPYDESCRIPTIONREFERENCE.ABSTRACTTEXT,
                    THERAPYDESCRIPTIONREFERENCE.URL,
                    THERAPYDESCRIPTIONREFERENCE.JOURNAL,
                    THERAPYDESCRIPTIONREFERENCE.AUTHORS,
                    THERAPYDESCRIPTIONREFERENCE.VOLUME,
                    THERAPYDESCRIPTIONREFERENCE.ISSUE,
                    THERAPYDESCRIPTIONREFERENCE.DATE,
                    THERAPYDESCRIPTIONREFERENCE.YEAR)
                    .values(id,
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
}

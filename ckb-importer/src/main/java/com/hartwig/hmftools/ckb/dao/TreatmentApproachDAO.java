package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.Tables.TREATMENTAPPROACHREFERENCE;
import static com.hartwig.hmftools.ckb.database.Tables.TREATMENTAPPROACHDRUGCLASS;

import static com.hartwig.hmftools.ckb.database.tables.Treatmentapproach.TREATMENTAPPROACH;

import com.hartwig.hmftools.ckb.datamodel.drug.DrugClass;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproaches.RelevantTreatmentApproaches;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;

class TreatmentApproachDAO {

    @NotNull
    private final DSLContext context;

    public TreatmentApproachDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAll() {
        // Note that deletions should go from branch to root
        context.deleteFrom(TREATMENTAPPROACHDRUGCLASS).execute();
        context.deleteFrom(TREATMENTAPPROACHREFERENCE).execute();

        context.deleteFrom(TREATMENTAPPROACH).execute();

    }

    public int write(@NotNull RelevantTreatmentApproaches treatmentApproaches) {
        int id = context.insertInto(TREATMENTAPPROACH,
                        TREATMENTAPPROACH.TREATMENTAPPROACHID,
                        TREATMENTAPPROACH.CREATEDATE,
                        TREATMENTAPPROACH.UPDATEDATE)
                .values(treatmentApproaches.id(), treatmentApproaches.createDate(), treatmentApproaches.updateDate())
                .returning(TREATMENTAPPROACH.ID)
                .fetchOne()
                .getValue(TREATMENTAPPROACH.ID);

        if (treatmentApproaches.drugClass() != null) {
            writeTreatmentDrugClass(treatmentApproaches.drugClass(), id);
        }

        for (Reference reference : treatmentApproaches.references()) {
            writeTreatmentReference(reference, id);
        }

        return id;
    }

    private void writeTreatmentDrugClass(@Nullable DrugClass drugClass, int treatmentApproachDrugClassId) {
        context.insertInto(TREATMENTAPPROACHDRUGCLASS,
                        TREATMENTAPPROACHDRUGCLASS.TREATMENTAPPROACHID,
                        TREATMENTAPPROACHDRUGCLASS.DRUGCLASSID,
                        TREATMENTAPPROACHDRUGCLASS.CREATEDATE,
                        TREATMENTAPPROACHDRUGCLASS.DRUGCLASS)
                .values(treatmentApproachDrugClassId, drugClass.id(), drugClass.createDate(), drugClass.drugClass())
                .execute();
    }

    private void writeTreatmentReference(@NotNull Reference reference, int treatmentApproachReferenceId) {
        context.insertInto(TREATMENTAPPROACHREFERENCE,
                        TREATMENTAPPROACHREFERENCE.TREATMENTAPPROACHID,
                        TREATMENTAPPROACHREFERENCE.REFERENCEID,
                        TREATMENTAPPROACHREFERENCE.PUBMEDID,
                        TREATMENTAPPROACHREFERENCE.TITLE,
                        TREATMENTAPPROACHREFERENCE.ABSTRACTTEXT,
                        TREATMENTAPPROACHREFERENCE.URL,
                        TREATMENTAPPROACHREFERENCE.JOURNAL,
                        TREATMENTAPPROACHREFERENCE.AUTHORS,
                        TREATMENTAPPROACHREFERENCE.VOLUME,
                        TREATMENTAPPROACHREFERENCE.ISSUE,
                        TREATMENTAPPROACHREFERENCE.DATE,
                        TREATMENTAPPROACHREFERENCE.YEAR)
                .values(treatmentApproachReferenceId,
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

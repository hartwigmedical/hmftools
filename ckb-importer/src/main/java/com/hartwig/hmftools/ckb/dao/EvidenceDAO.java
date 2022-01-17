package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.tables.Evidence.EVIDENCE;
import static com.hartwig.hmftools.ckb.database.tables.Evidencereference.EVIDENCEREFERENCE;
import static com.hartwig.hmftools.ckb.database.tables.Indicationevidence.INDICATIONEVIDENCE;
import static com.hartwig.hmftools.ckb.database.tables.Therapyevidence.THERAPYEVIDENCE;
import static com.hartwig.hmftools.ckb.database.tables.Treatmentapproachevidence.TREATMENTAPPROACHEVIDENCE;

import com.hartwig.hmftools.ckb.datamodel.evidence.Evidence;
import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproaches.RelevantTreatmentApproaches;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class EvidenceDAO {

    @NotNull
    private final DSLContext context;
    @NotNull
    private final TherapyDAO therapyDAO;
    @NotNull
    private final IndicationDAO indicationDAO;
    @NotNull
    private final TreatmentApproachDAO treatmentApproachDAO;

    public EvidenceDAO(@NotNull final DSLContext context, @NotNull final TherapyDAO therapyDAO,
            @NotNull final IndicationDAO indicationDAO, @NotNull final TreatmentApproachDAO treatmentApproachDAO) {
        this.context = context;
        this.therapyDAO = therapyDAO;
        this.indicationDAO = indicationDAO;
        this.treatmentApproachDAO = treatmentApproachDAO;
    }

    public void deleteAll() {
        // Note that deletions should go from branch to root
        context.deleteFrom(EVIDENCEREFERENCE).execute();
        context.deleteFrom(THERAPYEVIDENCE).execute();
        context.deleteFrom(INDICATIONEVIDENCE).execute();
        context.deleteFrom(TREATMENTAPPROACHEVIDENCE).execute();
        context.deleteFrom(EVIDENCE).execute();
    }

    public void write(@NotNull Evidence evidence, int ckbEntryId) {
        int id = context.insertInto(EVIDENCE,
                        EVIDENCE.CKBENTRYID,
                        EVIDENCE.CKBEVIDENCEID,
                        EVIDENCE.RESPONSETYPE,
                        EVIDENCE.EVIDENCETYPE,
                        EVIDENCE.EFFICACYEVIDENCE,
                        EVIDENCE.APPROVALSTATUS,
                        EVIDENCE.AMPCAPASCOEVIDENCELEVEL,
                        EVIDENCE.AMPCAPASCOINFERREDTIER)
                .values(ckbEntryId,
                        evidence.id(),
                        evidence.responseType(),
                        evidence.evidenceType(),
                        evidence.efficacyEvidence(),
                        evidence.approvalStatus(),
                        evidence.ampCapAscoEvidenceLevel(),
                        evidence.ampCapAscoInferredTier())
                .returning(EVIDENCE.ID)
                .fetchOne()
                .getValue(EVIDENCE.ID);

        int therapyId = therapyDAO.write(evidence.therapy());
        context.insertInto(THERAPYEVIDENCE, THERAPYEVIDENCE.EVIDENCEID, THERAPYEVIDENCE.THERAPYID).values(id, therapyId).execute();

        int indicationId = indicationDAO.write(evidence.indication());
        context.insertInto(INDICATIONEVIDENCE, INDICATIONEVIDENCE.EVIDENCEID, INDICATIONEVIDENCE.INDICATIONID)
                .values(id, indicationId)
                .execute();

        for (RelevantTreatmentApproaches treatmentApproaches : evidence.relevantTreatmentApproaches()) {
            int treatmentApproachId = treatmentApproachDAO.write(treatmentApproaches);

            context.insertInto(TREATMENTAPPROACHEVIDENCE,
                    TREATMENTAPPROACHEVIDENCE.EVIDENCEID,
                    TREATMENTAPPROACHEVIDENCE.TREATMENTAPPROACHEVIDENCEID).values(id, treatmentApproachId).execute();
        }

        for (Reference reference : evidence.references()) {
            writeReference(reference, id);
        }
    }

    private void writeReference(@NotNull Reference reference, int evidenceId) {
        context.insertInto(EVIDENCEREFERENCE,
                        EVIDENCEREFERENCE.EVIDENCEID,
                        EVIDENCEREFERENCE.CKBREFERENCEID,
                        EVIDENCEREFERENCE.PUBMEDID,
                        EVIDENCEREFERENCE.TITLE,
                        EVIDENCEREFERENCE.ABSTRACTTEXT,
                        EVIDENCEREFERENCE.URL,
                        EVIDENCEREFERENCE.JOURNAL,
                        EVIDENCEREFERENCE.AUTHORS,
                        EVIDENCEREFERENCE.VOLUME,
                        EVIDENCEREFERENCE.ISSUE,
                        EVIDENCEREFERENCE.DATE,
                        EVIDENCEREFERENCE.YEAR)
                .values(evidenceId,
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

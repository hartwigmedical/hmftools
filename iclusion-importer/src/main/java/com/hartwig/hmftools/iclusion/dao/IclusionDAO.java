package com.hartwig.hmftools.iclusion.dao;

import java.util.List;

import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;

import static com.hartwig.hmftools.iclusion.database.tables.Study.STUDY;
import static com.hartwig.hmftools.iclusion.database.tables.Studytumorlocations.STUDYTUMORLOCATIONS;
import static com.hartwig.hmftools.iclusion.database.tables.Studyblacklistedtumorlocations.STUDYBLACKLISTEDTUMORLOCATIONS;
import static com.hartwig.hmftools.iclusion.database.tables.Studymutationconditions.STUDYMUTATIONCONDITIONS;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class IclusionDAO {

    @NotNull
    private final DSLContext context;

    IclusionDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void write(@NotNull List<IclusionTrial> trials) {

        for (IclusionTrial trial : trials) {
            int id = context.insertInto(STUDY, STUDY.IDDB, STUDY.ACRONYM, STUDY.TITLE, STUDY.EUDRA, STUDY.NCT, STUDY.IPN, STUDY.CCMO)
                    .values(trial.id(), trial.acronym(), trial.title(), trial.eudra(), trial.nct(), trial.ipn(), trial.ccmo())
                    .returning(STUDY.ID)
                    .fetchOne()
                    .getValue(STUDY.ID);

            for (IclusionTumorLocation tumorLocation : trial.tumorLocations()) {
                context.insertInto(STUDYTUMORLOCATIONS, STUDYTUMORLOCATIONS.TUMORLOCATIONID, STUDYTUMORLOCATIONS.TUMORLOCATION)
                        .values(id, tumorLocation.primaryTumorLocation())
                        .execute();
            }

            for (IclusionTumorLocation blacklistedTumorLocation : trial.blacklistedTumorLocations()) {
                context.insertInto(STUDYBLACKLISTEDTUMORLOCATIONS,
                                STUDYBLACKLISTEDTUMORLOCATIONS.BLACKLISTEDTUMORLOCATIONID,
                                STUDYBLACKLISTEDTUMORLOCATIONS.BLACKLISTEDTUMORLOCATION)
                        .values(id, blacklistedTumorLocation.primaryTumorLocation())
                        .execute();
            }

            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                for (IclusionMutation mutation : mutationCondition.mutations()) {
                    context.insertInto(STUDYMUTATIONCONDITIONS,
                            STUDYMUTATIONCONDITIONS.MUTATIONCONDITIONID,
                            STUDYMUTATIONCONDITIONS.GENE,
                            STUDYMUTATIONCONDITIONS.MUTATION).values(id, mutation.gene(), mutation.name()).execute();
                }
            }
        }
    }
}
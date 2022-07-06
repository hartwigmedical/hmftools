package com.hartwig.hmftools.iclusion.dao;

import java.util.List;

import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;

import static com.hartwig.hmftools.iclusion.database.tables.Study.STUDY;
import static com.hartwig.hmftools.iclusion.database.tables.Tumorlocation.TUMORLOCATION;
import static com.hartwig.hmftools.iclusion.database.tables.Blacklistedtumorlocation.BLACKLISTEDTUMORLOCATION;
import static com.hartwig.hmftools.iclusion.database.tables.Mutationcondition.MUTATIONCONDITION;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

public class IclusionDAO {
    private static final Logger LOGGER = LogManager.getLogger(IclusionDAO.class);

    @NotNull
    private final DSLContext context;

    IclusionDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAll() {
        LOGGER.info("Deleting all data from iClusion database");
        // Note that deletions should go from branch to root
        context.deleteFrom(MUTATIONCONDITION).execute();
        context.deleteFrom(BLACKLISTEDTUMORLOCATION).execute();
        context.deleteFrom(TUMORLOCATION).execute();
        context.deleteFrom(STUDY).execute();
    }

    void write(@NotNull List<IclusionTrial> trials) {

        deleteAll();

        for (IclusionTrial trial : trials) {
            int id = context.insertInto(STUDY, STUDY.IDDB, STUDY.ACRONYM, STUDY.TITLE, STUDY.EUDRA, STUDY.NCT, STUDY.IPN, STUDY.CCMO)
                    .values(trial.id(),
                            trial.acronym(),
                            trial.title(),
                            trial.eudra().isEmpty() ? null : trial.eudra(),
                            trial.nct().isEmpty() ? null : trial.nct(),
                            trial.ipn().isEmpty() ? null : trial.ipn(),
                            trial.ccmo().isEmpty() ? null :  trial.ccmo())
                    .returning(STUDY.ID)
                    .fetchOne()
                    .getValue(STUDY.ID);

            for (IclusionTumorLocation tumorLocation : trial.tumorLocations()) {
                context.insertInto(TUMORLOCATION, TUMORLOCATION.TUMORLOCATIONID, TUMORLOCATION.TUMORLOCATION_)
                        .values(id, tumorLocation.primaryTumorLocation())
                        .execute();
            }

            for (IclusionTumorLocation blacklistedTumorLocation : trial.blacklistedTumorLocations()) {
                context.insertInto(BLACKLISTEDTUMORLOCATION,
                                BLACKLISTEDTUMORLOCATION.BLACKLISTEDTUMORLOCATIONID,
                                BLACKLISTEDTUMORLOCATION.BLACKLISTEDTUMORLOCATION_)
                        .values(id, blacklistedTumorLocation.primaryTumorLocation())
                        .execute();
            }

            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                for (IclusionMutation mutation : mutationCondition.mutations()) {
                    context.insertInto(MUTATIONCONDITION,
                            MUTATIONCONDITION.MUTATIONCONDITIONID,
                            MUTATIONCONDITION.GENE,
                            MUTATIONCONDITION.MUTATION).values(id, mutation.gene(), mutation.name()).execute();
                }
            }
        }
    }
}
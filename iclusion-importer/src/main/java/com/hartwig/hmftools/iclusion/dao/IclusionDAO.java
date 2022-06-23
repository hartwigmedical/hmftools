package com.hartwig.hmftools.iclusion.dao;

import java.util.List;

import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;

import static com.hartwig.hmftools.iclusion.database.tables.Study.STUDY;
import static com.hartwig.hmftools.iclusion.database.tables.Tumorlocations.TUMORLOCATIONS;
import static com.hartwig.hmftools.iclusion.database.tables.Blacklistedtumorlocations.BLACKLISTEDTUMORLOCATIONS;
import static com.hartwig.hmftools.iclusion.database.tables.Mutationconditions.MUTATIONCONDITIONS;

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
        context.deleteFrom(MUTATIONCONDITIONS).execute();
        context.deleteFrom(BLACKLISTEDTUMORLOCATIONS).execute();
        context.deleteFrom(TUMORLOCATIONS).execute();
        context.deleteFrom(STUDY).execute();
    }

    void write(@NotNull List<IclusionTrial> trials) {

        deleteAll();

        for (IclusionTrial trial : trials) {
            int id = context.insertInto(STUDY, STUDY.IDDB, STUDY.ACRONYM, STUDY.TITLE, STUDY.EUDRA, STUDY.NCT, STUDY.IPN, STUDY.CCMO)
                    .values(trial.id(), trial.acronym(), trial.title(), trial.eudra(), trial.nct(), trial.ipn(), trial.ccmo())
                    .returning(STUDY.ID)
                    .fetchOne()
                    .getValue(STUDY.ID);

            for (IclusionTumorLocation tumorLocation : trial.tumorLocations()) {
                context.insertInto(TUMORLOCATIONS, TUMORLOCATIONS.TUMORLOCATIONID, TUMORLOCATIONS.TUMORLOCATION)
                        .values(id, tumorLocation.primaryTumorLocation())
                        .execute();
            }

            for (IclusionTumorLocation blacklistedTumorLocation : trial.blacklistedTumorLocations()) {
                context.insertInto(BLACKLISTEDTUMORLOCATIONS,
                                BLACKLISTEDTUMORLOCATIONS.BLACKLISTEDTUMORLOCATIONID,
                                BLACKLISTEDTUMORLOCATIONS.BLACKLISTEDTUMORLOCATION)
                        .values(id, blacklistedTumorLocation.primaryTumorLocation())
                        .execute();
            }

            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                for (IclusionMutation mutation : mutationCondition.mutations()) {
                    context.insertInto(MUTATIONCONDITIONS,
                            MUTATIONCONDITIONS.MUTATIONCONDITIONID,
                            MUTATIONCONDITIONS.GENE,
                            MUTATIONCONDITIONS.MUTATION).values(id, mutation.gene(), mutation.name()).execute();
                }
            }
        }
    }
}
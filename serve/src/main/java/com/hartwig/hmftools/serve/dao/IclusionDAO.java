package com.hartwig.hmftools.serve.dao;

import java.sql.Timestamp;
import java.util.Date;
import java.util.List;

import com.hartwig.hmftools.iclusion.datamodel.IclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.IclusionMutationCondition;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;

import static com.hartwig.hmftools.serve.database.tables.Iclusionstudy.ICLUSIONSTUDY;
import static com.hartwig.hmftools.serve.database.tables.Iclusionstudytumorlocations.ICLUSIONSTUDYTUMORLOCATIONS;
import static com.hartwig.hmftools.serve.database.tables.Iclusionstudyblacklistedtumorlocations.ICLUSIONSTUDYBLACKLISTEDTUMORLOCATIONS;
import static com.hartwig.hmftools.serve.database.tables.Iclusionstudymutationconditions.ICLUSIONSTUDYMUTATIONCONDITIONS;

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
            int id = context.insertInto(ICLUSIONSTUDY,
                            ICLUSIONSTUDY.IDDB,
                            ICLUSIONSTUDY.ACRONYM,
                            ICLUSIONSTUDY.TITLE,
                            ICLUSIONSTUDY.EUDRA,
                            ICLUSIONSTUDY.NCT,
                            ICLUSIONSTUDY.IPN,
                            ICLUSIONSTUDY.CCMO)
                    .values(trial.id(), trial.acronym(), trial.title(), trial.eudra(), trial.nct(), trial.ipn(), trial.ccmo())
                    .returning(ICLUSIONSTUDY.ID)
                    .fetchOne()
                    .getValue(ICLUSIONSTUDY.ID);

            for (IclusionTumorLocation tumorLocation : trial.tumorLocations()) {
                context.insertInto(ICLUSIONSTUDYTUMORLOCATIONS,
                        ICLUSIONSTUDYTUMORLOCATIONS.TUMORLOCATIONID,
                        ICLUSIONSTUDYTUMORLOCATIONS.TUMORLOCATION).values(id, tumorLocation.primaryTumorLocation()).execute();
            }

            for (IclusionTumorLocation blacklistedTumorLocation : trial.blacklistedTumorLocations()) {
                context.insertInto(ICLUSIONSTUDYBLACKLISTEDTUMORLOCATIONS,
                                ICLUSIONSTUDYBLACKLISTEDTUMORLOCATIONS.BLACKLISTEDTUMORLOCATIONID,
                                ICLUSIONSTUDYBLACKLISTEDTUMORLOCATIONS.BLACKLISTEDTUMORLOCATION)
                        .values(id, blacklistedTumorLocation.primaryTumorLocation())
                        .execute();
            }

            for (IclusionMutationCondition mutationCondition : trial.mutationConditions()) {
                for (IclusionMutation mutation : mutationCondition.mutations()) {
                    context.insertInto(ICLUSIONSTUDYMUTATIONCONDITIONS,
                            ICLUSIONSTUDYMUTATIONCONDITIONS.MUTATIONCONDITIONID,
                            ICLUSIONSTUDYMUTATIONCONDITIONS.GENE,
                            ICLUSIONSTUDYMUTATIONCONDITIONS.MUTATION).values(id, mutation.gene(), mutation.name()).execute();
                }
            }
        }
    }
}
package com.hartwig.hmftools.ckb.dao;

import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;

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
        // TODO
    }

    public void write(@NotNull ClinicalTrial clinicalTrial, int ckbEntryId) {

    }
}

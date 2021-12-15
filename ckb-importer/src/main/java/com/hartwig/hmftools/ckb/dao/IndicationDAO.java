package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.tables.Indication.INDICATION;
import static com.hartwig.hmftools.ckb.database.tables.Indicationaltid.INDICATIONALTID;

import com.hartwig.hmftools.ckb.datamodel.indication.Indication;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class IndicationDAO {

    @NotNull
    private final DSLContext context;

    public IndicationDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAll() {
        // Note that deletions should go from branch to root
        context.deleteFrom(INDICATIONALTID).execute();
        context.deleteFrom(INDICATION).execute();
    }

    public int write(@NotNull Indication indication) {
        int id = context.insertInto(INDICATION,
                INDICATION.CKBINDICATIONID,
                INDICATION.NAME,
                INDICATION.SOURCE,
                INDICATION.DEFINITION,
                INDICATION.CURRENTPREFERREDTERM,
                INDICATION.LASTUPDATEDATEFROMDO,
                INDICATION.TERMID)
                .values(indication.id(),
                        indication.name(),
                        indication.source(),
                        indication.definition(),
                        indication.currentPreferredTerm(),
                        indication.lastUpdateDateFromDO(),
                        indication.termId())
                .returning(INDICATION.ID)
                .fetchOne()
                .getValue(INDICATION.ID);

        for (String altId : indication.altIds()) {
            context.insertInto(INDICATIONALTID, INDICATIONALTID.INDICATIONID, INDICATIONALTID.ALTID).values(id, altId).execute();
        }

        return id;
    }
}

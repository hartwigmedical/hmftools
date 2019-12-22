package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSINDICATION;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSMOLECULARPROFILE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSTHERAPY;

import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsIndication;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsTherapy;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class JaxTrialsDAOFunctions {

    private JaxTrialsDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull JaxTrials jaxTrials) {
        int id = context.insertInto(JAXTRIALS,
                JAXTRIALS.NCTID,
                JAXTRIALS.TITLE,
                JAXTRIALS.VARIANTREQUIREMENTS,
                JAXTRIALS.GENDER,
                JAXTRIALS.RECRUITMENT,
                JAXTRIALS.PHASE,
                JAXTRIALS.SPONSORS,
                JAXTRIALS.UPDATEDATE,
                JAXTRIALS.VICCENTRYID)
                .values(jaxTrials.nctId(),
                        jaxTrials.title(),
                        jaxTrials.variantRequirements(),
                        jaxTrials.gender(),
                        jaxTrials.recruitment(),
                        jaxTrials.phase(),
                        jaxTrials.sponsors(),
                        jaxTrials.updateDate(),
                        viccEntryId)
                .returning(JAXTRIALS.ID)
                .fetchOne()
                .getValue(JAXTRIALS.ID);

        for (JaxTrialsMolecularProfile molecularProfile : jaxTrials.molecularProfiles()) {
            context.insertInto(JAXTRIALSMOLECULARPROFILE,
                    JAXTRIALSMOLECULARPROFILE.REQUIREMENTTYPE,
                    JAXTRIALSMOLECULARPROFILE.PROFILENAME,
                    JAXTRIALSMOLECULARPROFILE.IDMOLECULARPROFILE,
                    JAXTRIALSMOLECULARPROFILE.JAXTRIALSID)
                    .values(molecularProfile.requirementType(), molecularProfile.profileName(), molecularProfile.id(), id)
                    .execute();
        }

        for (JaxTrialsIndication indication : jaxTrials.indications()) {
            context.insertInto(JAXTRIALSINDICATION,
                    JAXTRIALSINDICATION.NAME,
                    JAXTRIALSINDICATION.SOURCE,
                    JAXTRIALSINDICATION.IDINDICATION,
                    JAXTRIALSINDICATION.JAXTRIALSID).values(indication.name(), indication.source(), indication.id(), id).execute();
        }

        for (JaxTrialsTherapy therapy : jaxTrials.therapies()) {
            context.insertInto(JAXTRIALSTHERAPY, JAXTRIALSTHERAPY.THERAPYNAME, JAXTRIALSTHERAPY.IDTHERAPY, JAXTRIALSTHERAPY.JAXTRIALSID)
                    .values(therapy.therapyName(), therapy.id(), id)
                    .execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(JAXTRIALSMOLECULARPROFILE).execute();
        context.deleteFrom(JAXTRIALSINDICATION).execute();
        context.deleteFrom(JAXTRIALSTHERAPY).execute();

        context.deleteFrom(JAXTRIALS).execute();
    }
}

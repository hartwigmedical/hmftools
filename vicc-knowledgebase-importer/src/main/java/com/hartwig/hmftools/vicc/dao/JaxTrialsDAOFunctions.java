package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSINDICATION;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSMOLECULARPROFILE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSTHERAPY;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSVARIANTREQUIREMENTDETAILS;

import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsIndication;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsTherapy;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsVariantRequirementDetails;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class JaxTrialsDAOFunctions {

    private JaxTrialsDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull JaxTrials jaxTrials) {
        int id = context.insertInto(JAXTRIALS,
                JAXTRIALS.TITLE,
                JAXTRIALS.GENDER,
                JAXTRIALS.NCTID,
                JAXTRIALS.SPONSORS,
                JAXTRIALS.RECRUITMENT,
                JAXTRIALS.VARIANTREQUIREMENTS,
                JAXTRIALS.UPDATEDATE,
                JAXTRIALS.PHASE,
                JAXTRIALS.VICCENTRYID)
                .values(jaxTrials.title(),
                        jaxTrials.gender(),
                        jaxTrials.nctId(),
                        jaxTrials.sponsors(),
                        jaxTrials.recruitment(),
                        jaxTrials.variantRequirements(),
                        jaxTrials.updateDate(),
                        jaxTrials.updateDate(),
                        viccEntryId)
                .returning(JAXTRIALS.ID)
                .fetchOne()
                .getValue(JAXTRIALS.ID);

        for (JaxTrialsIndication indication : jaxTrials.indications()) {
            context.insertInto(JAXTRIALSINDICATION,
                    JAXTRIALSINDICATION.SOURCE,
                    JAXTRIALSINDICATION.IDINDICATION,
                    JAXTRIALSINDICATION.NAME,
                    JAXTRIALSINDICATION.JAXTRIALSID).values(indication.source(), indication.id(), indication.name(), id).execute();
        }

        for (JaxTrialsVariantRequirementDetails variantRequirementDetails : jaxTrials.variantRequirementDetails()) {
            int id1 = context.insertInto(JAXTRIALSVARIANTREQUIREMENTDETAILS,
                    JAXTRIALSVARIANTREQUIREMENTDETAILS.REQUIREMENTTYPE,
                    JAXTRIALSVARIANTREQUIREMENTDETAILS.JAXTRIALSID)
                    .values(variantRequirementDetails.requirementType(), id)
                    .returning(JAXTRIALSVARIANTREQUIREMENTDETAILS.ID)
                    .fetchOne()
                    .getValue(JAXTRIALSVARIANTREQUIREMENTDETAILS.ID);

            for (JaxTrialsMolecularProfile molecularProfile : variantRequirementDetails.molecularProfiles()) {
                context.insertInto(JAXTRIALSMOLECULARPROFILE,
                        JAXTRIALSMOLECULARPROFILE.PROFILENAME,
                        JAXTRIALSMOLECULARPROFILE.IDMOLECULARPROFILE,
                        JAXTRIALSMOLECULARPROFILE.JAXTRIALSVARIANTREQUIREMENTDETAILSID)
                        .values(molecularProfile.profileName(), molecularProfile.id(), id1)
                        .execute();
            }
        }

        for (JaxTrialsTherapy therapy : jaxTrials.therapies()) {
            context.insertInto(JAXTRIALSTHERAPY,
                    JAXTRIALSTHERAPY.IDTHERAPY,
                    JAXTRIALSTHERAPY.THERAPYNAME,
                    JAXTRIALSTHERAPY.JAXTRIALSID).values(therapy.id(), therapy.therapyName(), id).execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(JAXTRIALSMOLECULARPROFILE).execute();

        context.deleteFrom(JAXTRIALSINDICATION).execute();
        context.deleteFrom(JAXTRIALSVARIANTREQUIREMENTDETAILS).execute();
        context.deleteFrom(JAXTRIALSTHERAPY).execute();

        context.deleteFrom(JAXTRIALS).execute();
    }
}

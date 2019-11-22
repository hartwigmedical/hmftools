package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSINDICATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSMOLECULARPROFILE;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSTHERAPIES;
import static com.hartwig.hmftools.vicc.database.Tables.JAXTRIALSVARIANTREQUIREMENTDETAILS;

import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsIndications;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.jaxtrials.JaxTrialsTherapies;
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

        for (JaxTrialsIndications indications : jaxTrials.indications()) {
            context.insertInto(JAXTRIALSINDICATIONS,
                    JAXTRIALSINDICATIONS.SOURCE,
                    JAXTRIALSINDICATIONS.IDINDICATIONS,
                    JAXTRIALSINDICATIONS.NAME,
                    JAXTRIALSINDICATIONS.JAXTRIALSID).values(indications.source(), indications.id(), indications.name(), id).execute();
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

        for (JaxTrialsTherapies therapies : jaxTrials.therapies()) {
            context.insertInto(JAXTRIALSTHERAPIES,
                    JAXTRIALSTHERAPIES.IDTHERAPIES,
                    JAXTRIALSTHERAPIES.THERAPYNAME,
                    JAXTRIALSTHERAPIES.JAXTRIALSID).values(therapies.id(), therapies.therapyName(), id).execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(JAXTRIALSMOLECULARPROFILE).execute();

        context.deleteFrom(JAXTRIALSINDICATIONS).execute();
        context.deleteFrom(JAXTRIALSVARIANTREQUIREMENTDETAILS).execute();
        context.deleteFrom(JAXTRIALSTHERAPIES).execute();

        context.deleteFrom(JAXTRIALS).execute();
    }
}

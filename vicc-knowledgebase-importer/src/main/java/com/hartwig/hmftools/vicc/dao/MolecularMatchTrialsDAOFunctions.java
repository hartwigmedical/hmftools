package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSALTERATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSCONTACT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSGEO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSINTERVATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSLOCATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSLOCATIONS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOTHERGROUPLABEL;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOTHERNAME;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOVERALLCONTACT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSTAGS;

import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsIntervation;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsLocations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsTags;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class MolecularMatchTrialsDAOFunctions {

    private MolecularMatchTrialsDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull MolecularMatchTrials molecularMatchTrials) {
        int id = context.insertInto(MOLECULARMATCHTRIALS,
                MOLECULARMATCHTRIALS.STATUS,
                MOLECULARMATCHTRIALS.STARTDATE,
                MOLECULARMATCHTRIALS.TITLE,
                MOLECULARMATCHTRIALS.SCORE,
                MOLECULARMATCHTRIALS.BRIEFTITLE,
                MOLECULARMATCHTRIALS.LINK,
                MOLECULARMATCHTRIALS.PHASE,
                MOLECULARMATCHTRIALS.IDMOLECULARMATCHTRIALS,
                MOLECULARMATCHTRIALS.STUDYTYPE,
                MOLECULARMATCHTRIALS.VICCENTRYID)
                .values(molecularMatchTrials.status(),
                        molecularMatchTrials.startDate(),
                        molecularMatchTrials.title(),
                        molecularMatchTrials.score(),
                        molecularMatchTrials.briefTitle(),
                        molecularMatchTrials.link(),
                        molecularMatchTrials.phase(),
                        molecularMatchTrials.id(),
                        molecularMatchTrials.studyType(),
                        viccEntryId)
                .returning(MOLECULARMATCHTRIALS.ID)
                .fetchOne()
                .getValue(MOLECULARMATCHTRIALS.ID);

        for (String molecularAlterations : molecularMatchTrials.molecularAlterations()) {
            context.insertInto(MOLECULARMATCHTRIALSALTERATIONS,
                    MOLECULARMATCHTRIALSALTERATIONS.MOLECULARALTERATIONS,
                    MOLECULARMATCHTRIALSALTERATIONS.MOLECULARMATCHTRIALSID).values(molecularAlterations, id).execute();
        }

        for (MolecularMatchTrialsIntervation intervation : molecularMatchTrials.intervation()) {
            int idIntervation = context.insertInto(MOLECULARMATCHTRIALSINTERVATIONS,
                    MOLECULARMATCHTRIALSINTERVATIONS.INTERVENTION_NAME,
                    MOLECULARMATCHTRIALSINTERVATIONS.DESCRIPTION,
                    MOLECULARMATCHTRIALSINTERVATIONS.INTERVENTION_TYPE,
                    MOLECULARMATCHTRIALSINTERVATIONS.MOLECULARMATCHTRIALSID)
                    .values(intervation.intervention_name(), intervation.description(), intervation.intervention_type(), id)
                    .returning(MOLECULARMATCHTRIALSINTERVATIONS.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHTRIALSINTERVATIONS.ID);

            if (intervation.other_name() != null) {
                for (String otherName : intervation.other_name()) {
                    context.insertInto(MOLECULARMATCHTRIALSOTHERNAME,
                            MOLECULARMATCHTRIALSOTHERNAME.OTHER_NAME,
                            MOLECULARMATCHTRIALSOTHERNAME.MOLECULARMATCHTRIALSINTERVATIONSID).values(otherName, idIntervation).execute();
                }
            }

            if (intervation.arm_group_label() != null) {
                for (String armGroupLabel : intervation.arm_group_label()) {
                    context.insertInto(MOLECULARMATCHTRIALSOTHERGROUPLABEL,
                            MOLECULARMATCHTRIALSOTHERGROUPLABEL.ARM_GROUP_LABEL,
                            MOLECULARMATCHTRIALSOTHERGROUPLABEL.MOLECULARMATCHTRIALSINTERVATIONSID)
                            .values(armGroupLabel, idIntervation)
                            .execute();
                }
            }
        }

        for (MolecularMatchTrialsLocations locations : molecularMatchTrials.locations()) {
            int idLocations = context.insertInto(MOLECULARMATCHTRIALSLOCATIONS,
                    MOLECULARMATCHTRIALSLOCATIONS.STATUS,
                    MOLECULARMATCHTRIALSLOCATIONS.LAST_NAME,
                    MOLECULARMATCHTRIALSLOCATIONS.EMAIL,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.EMAIL_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.LAST_NAME_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE_EXT_BACKUP,
                    MOLECULARMATCHTRIALSLOCATIONS.PHONE_EXT,
                    MOLECULARMATCHTRIALSLOCATIONS.CITY,
                    MOLECULARMATCHTRIALSLOCATIONS.VALID,
                    MOLECULARMATCHTRIALSLOCATIONS.ZIP,
                    MOLECULARMATCHTRIALSLOCATIONS.CREATED,
                    MOLECULARMATCHTRIALSLOCATIONS.COUNTRY,
                    MOLECULARMATCHTRIALSLOCATIONS.NUMBER,
                    MOLECULARMATCHTRIALSLOCATIONS.IDLOCATIONS,
                    MOLECULARMATCHTRIALSLOCATIONS.LASTUPDATED,
                    MOLECULARMATCHTRIALSLOCATIONS.STATE,
                    MOLECULARMATCHTRIALSLOCATIONS.STREET,
                    MOLECULARMATCHTRIALSLOCATIONS.PO_BOX,
                    MOLECULARMATCHTRIALSLOCATIONS.FAILEDGEOCODE,
                    MOLECULARMATCHTRIALSLOCATIONS.VALIDMESSAGE,
                    MOLECULARMATCHTRIALSLOCATIONS.NAME,
                    MOLECULARMATCHTRIALSLOCATIONS.MOLECULARMATCHTRIALSID)
                    .values(locations.status(),
                            locations.last_name(),
                            locations.email(),
                            locations.phone(),
                            locations.phone_backup(),
                            locations.email_backup(),
                            locations.last_name_backup(),
                            locations.phone_ext_backup(),
                            locations.phone_ext(),
                            locations.city(),
                            locations.valid(),
                            locations.zip(),
                            locations.created(),
                            locations.country(),
                            locations.number(),
                            locations.id(),
                            locations.lastUpdated(),
                            locations.state(),
                            locations.street(),
                            locations.po_box(),
                            locations.failedGeocode(),
                            locations.validMessage(),
                            locations.name(),
                            id)
                    .returning(MOLECULARMATCHTRIALSLOCATIONS.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHTRIALSLOCATIONS.ID);

            if (locations.contact() != null) {
                context.insertInto(MOLECULARMATCHTRIALSCONTACT,
                        MOLECULARMATCHTRIALSCONTACT.PHONE,
                        MOLECULARMATCHTRIALSCONTACT.NAME,
                        MOLECULARMATCHTRIALSCONTACT.EMAIL,
                        MOLECULARMATCHTRIALSCONTACT.MOLECULARMATCHTRIALSLOCATIONSID)
                        .values(locations.contact().phone(), locations.contact().name(), locations.contact().email(), idLocations)
                        .execute();
            }

            if (locations.location() != null) {
                int idLocation = context.insertInto(MOLECULARMATCHTRIALSLOCATION,
                        MOLECULARMATCHTRIALSLOCATION.TYPE,
                        MOLECULARMATCHTRIALSLOCATION.MOLECULARMATCHTRIALSLOCATIONSID)
                        .values(locations.location().type(), idLocations)
                        .returning(MOLECULARMATCHTRIALSLOCATION.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHTRIALSLOCATION.ID);

                for (String coordinates : locations.location().coordinates()) {
                    context.insertInto(MOLECULARMATCHTRIALSCOORDINATES,
                            MOLECULARMATCHTRIALSCOORDINATES.COORDINATES,
                            MOLECULARMATCHTRIALSCOORDINATES.MOLECULARMATCHTRIALSLOCATIONID).values(coordinates, idLocation).execute();
                }
            }

            if (locations.geo() != null) {
                context.insertInto(MOLECULARMATCHTRIALSGEO,
                        MOLECULARMATCHTRIALSGEO.LAT,
                        MOLECULARMATCHTRIALSGEO.LON,
                        MOLECULARMATCHTRIALSGEO.MOLECULARMATCHTRIALSLOCATIONSID)
                        .values(locations.geo().lat(), locations.geo().lon(), idLocations)
                        .execute();
            }
        }

        if (molecularMatchTrials.overallContact() != null) {
            context.insertInto(MOLECULARMATCHTRIALSOVERALLCONTACT,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.PHONE,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.LAST_NAME,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.PHONE_EXT,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.COUNTRY,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.EMAIL,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.AFFILIATION,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.CITY,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.NAME,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.ZIP,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.URL,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.STREET,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.TYPE,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.MOLECULARMATCHTRIALSID)
                    .values(molecularMatchTrials.overallContact().phone(),
                            molecularMatchTrials.overallContact().last_name(),
                            molecularMatchTrials.overallContact().phone_ext(),
                            molecularMatchTrials.overallContact().country(),
                            molecularMatchTrials.overallContact().email(),
                            molecularMatchTrials.overallContact().affiliation(),
                            molecularMatchTrials.overallContact().city(),
                            molecularMatchTrials.overallContact().name(),
                            molecularMatchTrials.overallContact().zip(),
                            molecularMatchTrials.overallContact().url(),
                            molecularMatchTrials.overallContact().street(),
                            molecularMatchTrials.overallContact().type(),
                            id)
                    .execute();
        }

        for (MolecularMatchTrialsTags tags : molecularMatchTrials.tags()) {
            context.insertInto(MOLECULARMATCHTRIALSTAGS,
                    MOLECULARMATCHTRIALSTAGS.FACET,
                    MOLECULARMATCHTRIALSTAGS.COMPOSITEKEY,
                    MOLECULARMATCHTRIALSTAGS.SUPPRESS,
                    MOLECULARMATCHTRIALSTAGS.FILTERTYPE,
                    MOLECULARMATCHTRIALSTAGS.TERM,
                    MOLECULARMATCHTRIALSTAGS.CUSTOM,
                    MOLECULARMATCHTRIALSTAGS.PRIORITY,
                    MOLECULARMATCHTRIALSTAGS.ALIAS,
                    MOLECULARMATCHTRIALSTAGS.MANUALSUPPRESS,
                    MOLECULARMATCHTRIALSTAGS.GENERATEDBY,
                    MOLECULARMATCHTRIALSTAGS.GENERATEDBYTERM,
                    MOLECULARMATCHTRIALSTAGS.IDTAGS,
                    MOLECULARMATCHTRIALSTAGS.MANUALPRIORITY,
                    MOLECULARMATCHTRIALSTAGS.MOLECULARMATCHTRIALSID)
                    .values(tags.facet(),
                            tags.compositeKey(),
                            tags.suppress(),
                            tags.filterType(),
                            tags.term(),
                            tags.custom(),
                            tags.priority(),
                            tags.alias(),
                            tags.manualSuppress(),
                            tags.generatedBy(),
                            tags.generatedByTerm(),
                            tags.id(),
                            tags.manualPriority(),
                            id)
                    .execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(MOLECULARMATCHTRIALS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSALTERATIONS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSINTERVATIONS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOTHERNAME).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOTHERGROUPLABEL).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOVERALLCONTACT).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSTAGS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSLOCATIONS).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSCONTACT).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSGEO).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSLOCATION).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSCOORDINATES).execute();
    }
}

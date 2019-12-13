package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALS;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSALTERATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSARMGROUPLABEL;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSCONTACT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSGEO;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSINTERVENTION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSLOCATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOTHERNAME;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSOVERALLCONTACT;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSSUBLOCATION;
import static com.hartwig.hmftools.vicc.database.Tables.MOLECULARMATCHTRIALSTAG;

import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsGeo;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsIntervention;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsOverallContact;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsSubLocation;
import com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials.MolecularMatchTrialsTag;

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
                MOLECULARMATCHTRIALS.BRIEFTITLE,
                MOLECULARMATCHTRIALS.STUDYTYPE,
                MOLECULARMATCHTRIALS.SCORE,
                MOLECULARMATCHTRIALS.LINK,
                MOLECULARMATCHTRIALS.PHASE,
                MOLECULARMATCHTRIALS.IDMOLECULARMATCHTRIALS,
                MOLECULARMATCHTRIALS.VICCENTRYID)
                .values(molecularMatchTrials.status(),
                        molecularMatchTrials.startDate(),
                        molecularMatchTrials.title(),
                        molecularMatchTrials.briefTitle(),
                        molecularMatchTrials.studyType(),
                        molecularMatchTrials.score(),
                        molecularMatchTrials.link(),
                        molecularMatchTrials.phase(),
                        molecularMatchTrials.id(),
                        viccEntryId)
                .returning(MOLECULARMATCHTRIALS.ID)
                .fetchOne()
                .getValue(MOLECULARMATCHTRIALS.ID);

        for (String molecularAlteration : molecularMatchTrials.molecularAlterations()) {
            context.insertInto(MOLECULARMATCHTRIALSALTERATION,
                    MOLECULARMATCHTRIALSALTERATION.MOLECULARALTERATION,
                    MOLECULARMATCHTRIALSALTERATION.MOLECULARMATCHTRIALSID).values(molecularAlteration, id).execute();
        }

        for (MolecularMatchTrialsIntervention intervention : molecularMatchTrials.interventions()) {
            int idIntervention = context.insertInto(MOLECULARMATCHTRIALSINTERVENTION,
                    MOLECULARMATCHTRIALSINTERVENTION.INTERVENTIONNAME,
                    MOLECULARMATCHTRIALSINTERVENTION.INTERVENTIONTYPE,
                    MOLECULARMATCHTRIALSINTERVENTION.DESCRIPTION,
                    MOLECULARMATCHTRIALSINTERVENTION.MOLECULARMATCHTRIALSID)
                    .values(intervention.interventionName(), intervention.interventionType(), intervention.description(), id)
                    .returning(MOLECULARMATCHTRIALSINTERVENTION.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHTRIALSINTERVENTION.ID);

            for (String otherName : intervention.otherNames()) {
                context.insertInto(MOLECULARMATCHTRIALSOTHERNAME,
                        MOLECULARMATCHTRIALSOTHERNAME.OTHERNAME,
                        MOLECULARMATCHTRIALSOTHERNAME.MOLECULARMATCHTRIALSINTERVENTIONID).values(otherName, idIntervention).execute();
            }

            for (String armGroupLabel : intervention.armGroupLabels()) {
                context.insertInto(MOLECULARMATCHTRIALSARMGROUPLABEL,
                        MOLECULARMATCHTRIALSARMGROUPLABEL.ARMGROUPLABEL,
                        MOLECULARMATCHTRIALSARMGROUPLABEL.MOLECULARMATCHTRIALSINTERVENTIONID)
                        .values(armGroupLabel, idIntervention)
                        .execute();
            }
        }

        for (MolecularMatchTrialsLocation location : molecularMatchTrials.locations()) {
            int idLocation = context.insertInto(MOLECULARMATCHTRIALSLOCATION,
                    MOLECULARMATCHTRIALSLOCATION.STATUS,
                    MOLECULARMATCHTRIALSLOCATION.NAME,
                    MOLECULARMATCHTRIALSLOCATION.LASTNAME,
                    MOLECULARMATCHTRIALSLOCATION.EMAIL,
                    MOLECULARMATCHTRIALSLOCATION.PHONE,
                    MOLECULARMATCHTRIALSLOCATION.PHONEEXT,
                    MOLECULARMATCHTRIALSLOCATION.LASTNAMEBACKUP,
                    MOLECULARMATCHTRIALSLOCATION.EMAILBACKUP,
                    MOLECULARMATCHTRIALSLOCATION.PHONEBACKUP,
                    MOLECULARMATCHTRIALSLOCATION.PHONEEXTBACKUP,
                    MOLECULARMATCHTRIALSLOCATION.STREET,
                    MOLECULARMATCHTRIALSLOCATION.ZIP,
                    MOLECULARMATCHTRIALSLOCATION.CITY,
                    MOLECULARMATCHTRIALSLOCATION.STATE,
                    MOLECULARMATCHTRIALSLOCATION.COUNTRY,
                    MOLECULARMATCHTRIALSLOCATION.NUMBER,
                    MOLECULARMATCHTRIALSLOCATION.POBOX,
                    MOLECULARMATCHTRIALSLOCATION.IDLOCATION,
                    MOLECULARMATCHTRIALSLOCATION.VALID,
                    MOLECULARMATCHTRIALSLOCATION.VALIDMESSAGE,
                    MOLECULARMATCHTRIALSLOCATION.CREATED,
                    MOLECULARMATCHTRIALSLOCATION.LASTUPDATED,
                    MOLECULARMATCHTRIALSLOCATION.FAILEDGEOCODE,
                    MOLECULARMATCHTRIALSLOCATION.MOLECULARMATCHTRIALSID)
                    .values(location.status(),
                            location.name(),
                            location.lastName(),
                            location.email(),
                            location.phone(),
                            location.phoneExt(),
                            location.lastNameBackup(),
                            location.emailBackup(),
                            location.phoneBackup(),
                            location.phoneExtBackup(),
                            location.street(),
                            location.zip(),
                            location.city(),
                            location.state(),
                            location.country(),
                            location.number(),
                            location.poBox(),
                            location.id(),
                            location.valid(),
                            location.validMessage(),
                            location.created(),
                            location.lastUpdated(),
                            location.failedGeocode(),
                            id)
                    .returning(MOLECULARMATCHTRIALSLOCATION.ID)
                    .fetchOne()
                    .getValue(MOLECULARMATCHTRIALSLOCATION.ID);

            MolecularMatchTrialsContact contact = location.contact();
            if (contact != null) {
                context.insertInto(MOLECULARMATCHTRIALSCONTACT,
                        MOLECULARMATCHTRIALSCONTACT.NAME,
                        MOLECULARMATCHTRIALSCONTACT.EMAIL,
                        MOLECULARMATCHTRIALSCONTACT.PHONE,
                        MOLECULARMATCHTRIALSCONTACT.MOLECULARMATCHTRIALSLOCATIONID)
                        .values(contact.name(), contact.email(), contact.phone(), idLocation)
                        .execute();
            }

            MolecularMatchTrialsSubLocation subLocation = location.subLocation();
            if (subLocation != null) {
                int idSubLocation = context.insertInto(MOLECULARMATCHTRIALSSUBLOCATION,
                        MOLECULARMATCHTRIALSSUBLOCATION.TYPE,
                        MOLECULARMATCHTRIALSSUBLOCATION.MOLECULARMATCHTRIALSLOCATIONID)
                        .values(subLocation.type(), idLocation)
                        .returning(MOLECULARMATCHTRIALSSUBLOCATION.ID)
                        .fetchOne()
                        .getValue(MOLECULARMATCHTRIALSSUBLOCATION.ID);

                for (String coordinate : subLocation.coordinates()) {
                    context.insertInto(MOLECULARMATCHTRIALSCOORDINATES,
                            MOLECULARMATCHTRIALSCOORDINATES.COORDINATES,
                            MOLECULARMATCHTRIALSCOORDINATES.MOLECULARMATCHTRIALSSUBLOCATIONID).values(coordinate, idSubLocation).execute();
                }
            }

            MolecularMatchTrialsGeo geo = location.geo();
            if (geo != null) {
                context.insertInto(MOLECULARMATCHTRIALSGEO,
                        MOLECULARMATCHTRIALSGEO.LAT,
                        MOLECULARMATCHTRIALSGEO.LON,
                        MOLECULARMATCHTRIALSGEO.MOLECULARMATCHTRIALSLOCATIONID)
                        .values(geo.lat(), geo.lon(), idLocation)
                        .execute();
            }
        }

        MolecularMatchTrialsOverallContact overallContact = molecularMatchTrials.overallContact();
        if (overallContact != null) {
            context.insertInto(MOLECULARMATCHTRIALSOVERALLCONTACT,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.NAME,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.TYPE,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.AFFILIATION,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.LASTNAME,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.EMAIL,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.PHONE,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.PHONEEXT,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.STREET,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.ZIP,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.CITY,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.COUNTRY,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.URL,
                    MOLECULARMATCHTRIALSOVERALLCONTACT.MOLECULARMATCHTRIALSID)
                    .values(overallContact.name(),
                            overallContact.type(),
                            overallContact.affiliation(),
                            overallContact.lastName(),
                            overallContact.email(),
                            overallContact.phone(),
                            overallContact.phoneExt(),
                            overallContact.street(),
                            overallContact.zip(),
                            overallContact.city(),
                            overallContact.country(),
                            overallContact.url(),
                            id)
                    .execute();
        }

        for (MolecularMatchTrialsTag tag : molecularMatchTrials.tags()) {
            context.insertInto(MOLECULARMATCHTRIALSTAG,
                    MOLECULARMATCHTRIALSTAG.FACET,
                    MOLECULARMATCHTRIALSTAG.COMPOSITEKEY,
                    MOLECULARMATCHTRIALSTAG.SUPPRESS,
                    MOLECULARMATCHTRIALSTAG.FILTERTYPE,
                    MOLECULARMATCHTRIALSTAG.TERM,
                    MOLECULARMATCHTRIALSTAG.CUSTOM,
                    MOLECULARMATCHTRIALSTAG.PRIORITY,
                    MOLECULARMATCHTRIALSTAG.ALIAS,
                    MOLECULARMATCHTRIALSTAG.MANUALSUPPRESS,
                    MOLECULARMATCHTRIALSTAG.GENERATEDBY,
                    MOLECULARMATCHTRIALSTAG.GENERATEDBYTERM,
                    MOLECULARMATCHTRIALSTAG.IDTAG,
                    MOLECULARMATCHTRIALSTAG.MANUALPRIORITY,
                    MOLECULARMATCHTRIALSTAG.MOLECULARMATCHTRIALSID)
                    .values(tag.facet(),
                            tag.compositeKey(),
                            tag.suppress(),
                            tag.filterType(),
                            tag.term(),
                            tag.custom(),
                            tag.priority(),
                            tag.alias(),
                            tag.manualSuppress(),
                            tag.generatedBy(),
                            tag.generatedByTerm(),
                            tag.id(),
                            tag.manualPriority(),
                            id)
                    .execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        // Tables that are part of a location
        context.deleteFrom(MOLECULARMATCHTRIALSCOORDINATES).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSSUBLOCATION).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSGEO).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSCONTACT).execute();

        // Tables that are part of an intervention
        context.deleteFrom(MOLECULARMATCHTRIALSOTHERNAME).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSARMGROUPLABEL).execute();

        // Tables that are part of the main entry object
        context.deleteFrom(MOLECULARMATCHTRIALSTAG).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSOVERALLCONTACT).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSLOCATION).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSINTERVENTION).execute();
        context.deleteFrom(MOLECULARMATCHTRIALSALTERATION).execute();

        context.deleteFrom(MOLECULARMATCHTRIALS).execute();
    }
}

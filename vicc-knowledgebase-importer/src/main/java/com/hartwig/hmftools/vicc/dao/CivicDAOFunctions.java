package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.CIVIC;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICASSERTIONS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCLINICALTRIAL;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCLINVARENTRIES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDESCRIPTION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDISEASE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDRUGS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICERROR;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSCLINICALTRIAL;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSPUBLICATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICHGVSEXPRESSIONS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDON;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONAVATARS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONORGANIZATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONPROFILEIMAGE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTCOMMENTEDONUSER;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIED;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDAVATARS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDORGANIZATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDPROFILEIMAGE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTMODIFIEDUSER;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWED;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDAVATARS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDORGANIZATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDPROFILEIMAGE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLASTREVIEWEDUSER;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICLIFECYCLEACTIONS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICPUBLICATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTALIASES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPSCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPSTYPES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTSGROUPSVARIANTS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTTYPES;

import com.hartwig.hmftools.vicc.datamodel.Civic;
import com.hartwig.hmftools.vicc.datamodel.CivicClinicalTrial;
import com.hartwig.hmftools.vicc.datamodel.CivicDrugs;
import com.hartwig.hmftools.vicc.datamodel.CivicEvidenceItems;
import com.hartwig.hmftools.vicc.datamodel.CivicSource;
import com.hartwig.hmftools.vicc.datamodel.CivicUser;
import com.hartwig.hmftools.vicc.datamodel.CivicVariantGroup;
import com.hartwig.hmftools.vicc.datamodel.CivicVariantTypes;
import com.hartwig.hmftools.vicc.datamodel.CivicVariants;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class CivicDAOFunctions {

    private CivicDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull Civic civic) {
        int id = context.insertInto(CIVIC,
                CIVIC.ENTREZNAME,
                CIVIC.CIVICACTIONABILITYSCORE,
                CIVIC.ALLELEREGISTRYID,
                CIVIC.GENEID,
                CIVIC.NAME,
                CIVIC.ENTREZID,
                CIVIC.TYPE,
                CIVIC.IDCIVIC,
                CIVIC.DESCRIPTION,
                CIVIC.VICCENTRYID)
                .values(civic.entrezName(),
                        civic.civicActionabilityScore(),
                        civic.alleleRegistryId(),
                        civic.geneId(),
                        civic.name(),
                        civic.entrezId(),
                        civic.type(),
                        civic.id(),
                        civic.description(),
                        viccEntryId)
                .returning(CIVIC.ID)
                .fetchOne()
                .getValue(CIVIC.ID);

        for (String assertions : civic.assertions()) {
            context.insertInto(CIVICASSERTIONS, CIVICASSERTIONS.ASSERTIONS, CIVICASSERTIONS.CIVICID).values(assertions, id).execute();
        }

        for (String hgvsExpression : civic.hgvs_expressions()) {
            context.insertInto(CIVICHGVSEXPRESSIONS, CIVICHGVSEXPRESSIONS.HGVS_EXPRESSIONS, CIVICHGVSEXPRESSIONS.CIVICID)
                    .values(hgvsExpression, id)
                    .execute();
        }

        for (String clinvarEntries : civic.clinvarEntries()) {
            context.insertInto(CIVICCLINVARENTRIES, CIVICCLINVARENTRIES.CLINVARENTRIES, CIVICCLINVARENTRIES.CIVICID)
                    .values(clinvarEntries, id)
                    .execute();
        }

        for (String variantAliases : civic.variantAliases()) {
            context.insertInto(CIVICVARIANTALIASES, CIVICVARIANTALIASES.VARIANTALIASES, CIVICVARIANTALIASES.CIVICID)
                    .values(variantAliases, id)
                    .execute();
        }

        for (CivicVariantTypes variantTypes : civic.variantTypes()) {
            context.insertInto(CIVICVARIANTTYPES,
                    CIVICVARIANTTYPES.DISPLAYNAME,
                    CIVICVARIANTTYPES.DESCRIPTION,
                    CIVICVARIANTTYPES.URL,
                    CIVICVARIANTTYPES.SOID,
                    CIVICVARIANTTYPES.IDVARIANTTYPES,
                    CIVICVARIANTTYPES.NAME,
                    CIVICVARIANTTYPES.CIVICID)
                    .values(variantTypes.displayName(),
                            variantTypes.description(),
                            variantTypes.url(),
                            variantTypes.soId(),
                            variantTypes.id(),
                            variantTypes.name(),
                            id)
                    .execute();
        }

        if (civic.provisional_values() != null) {
            context.insertInto(CIVICDESCRIPTION, CIVICDESCRIPTION.REVISIONID, CIVICDESCRIPTION.VALUE, CIVICDESCRIPTION.CIVICID)
                    .values(civic.provisional_values().revision_id(), civic.provisional_values().value(), id)
                    .execute();
        }

        context.insertInto(CIVICCOORDINATES,
                CIVICCOORDINATES.CHROMOSOME2,
                CIVICCOORDINATES.REFERENCEBASES,
                CIVICCOORDINATES.START2,
                CIVICCOORDINATES.VARIANTBASES,
                CIVICCOORDINATES.STOP,
                CIVICCOORDINATES.STOP2,
                CIVICCOORDINATES.REPRESENTATIVETRANSCRIPT2,
                CIVICCOORDINATES.START,
                CIVICCOORDINATES.REPRESENTATIVETRANSCRIPT,
                CIVICCOORDINATES.ENSEMBLVERSION,
                CIVICCOORDINATES.CHROMOSOME,
                CIVICCOORDINATES.REFERENCEBUILD,
                CIVICCOORDINATES.CIVICID)
                .values(civic.coordinates().chromosome2(),
                        civic.coordinates().referenceBases(),
                        civic.coordinates().start2(),
                        civic.coordinates().variantBases(),
                        civic.coordinates().stop(),
                        civic.coordinates().stop2(),
                        civic.coordinates().representativeTranscript2(),
                        civic.coordinates().start(),
                        civic.coordinates().representativeTranscript(),
                        civic.coordinates().ensemblVersion(),
                        civic.coordinates().chromosome(),
                        civic.coordinates().referenceBuild(),
                        id)
                .execute();

        for (CivicVariantGroup variantGroup : civic.variantGroups()) {
            int idVariantGroup = context.insertInto(CIVICVARIANTSGROUPS,
                    CIVICVARIANTSGROUPS.IDVARIANTSGROUPS,
                    CIVICVARIANTSGROUPS.TYPE,
                    CIVICVARIANTSGROUPS.DESCRIPTION,
                    CIVICVARIANTSGROUPS.NAME,
                    CIVICVARIANTSGROUPS.CIVICID)
                    .values(variantGroup.id(), variantGroup.type(), variantGroup.description(), variantGroup.name(), id)
                    .returning(CIVICVARIANTSGROUPS.ID)
                    .fetchOne()
                    .getValue(CIVICVARIANTSGROUPS.ID);

            for (CivicVariants variants : variantGroup.variants()) {
                int idVariantGroupVariants = context.insertInto(CIVICVARIANTSGROUPSVARIANTS,
                        CIVICVARIANTSGROUPSVARIANTS.ENTREZ_NAME,
                        CIVICVARIANTSGROUPSVARIANTS.DESCRIPTION,
                        CIVICVARIANTSGROUPSVARIANTS.CIVIC_ACTIONABILITY_SCORE,
                        CIVICVARIANTSGROUPSVARIANTS.GENE_ID,
                        CIVICVARIANTSGROUPSVARIANTS.ENTREZ_ID,
                        CIVICVARIANTSGROUPSVARIANTS.TYPE,
                        CIVICVARIANTSGROUPSVARIANTS.IDVARIANTS,
                        CIVICVARIANTSGROUPSVARIANTS.NAME,
                        CIVICVARIANTSGROUPSVARIANTS.CIVICVARIANTSGROUPSID)
                        .values(variants.entrez_name(),
                                variants.description(),
                                variants.civic_actionability_score(),
                                variants.gene_id(),
                                variants.entrez_id(),
                                variants.type(),
                                variants.id(),
                                variants.name(),
                                idVariantGroup)
                        .returning(CIVICVARIANTSGROUPSVARIANTS.ID)
                        .fetchOne()
                        .getValue(CIVICVARIANTSGROUPSVARIANTS.ID);

                if (variants.coordinates() != null) {
                    context.insertInto(CIVICVARIANTSGROUPSCOORDINATES,
                            CIVICVARIANTSGROUPSCOORDINATES.CHROMOSOME2,
                            CIVICVARIANTSGROUPSCOORDINATES.REFERENCEBASES,
                            CIVICVARIANTSGROUPSCOORDINATES.START2,
                            CIVICVARIANTSGROUPSCOORDINATES.VARIANTBASES,
                            CIVICVARIANTSGROUPSCOORDINATES.STOP,
                            CIVICVARIANTSGROUPSCOORDINATES.STOP2,
                            CIVICVARIANTSGROUPSCOORDINATES.REPRESENTATIVETRANSCRIPT2,
                            CIVICVARIANTSGROUPSCOORDINATES.START,
                            CIVICVARIANTSGROUPSCOORDINATES.REPRESENTATIVETRANSCRIPT,
                            CIVICVARIANTSGROUPSCOORDINATES.ENSEMBLVERSION,
                            CIVICVARIANTSGROUPSCOORDINATES.CHROMOSOME,
                            CIVICVARIANTSGROUPSCOORDINATES.REFERENCEBUILD,
                            CIVICVARIANTSGROUPSCOORDINATES.CIVICVARIANTSGROUPSVARIANTSID)
                            .values(variants.coordinates().chromosome2(),
                                    variants.coordinates().referenceBases(),
                                    variants.coordinates().start2(),
                                    variants.coordinates().variantBases(),
                                    variants.coordinates().stop(),
                                    variants.coordinates().stop2(),
                                    variants.coordinates().representativeTranscript2(),
                                    variants.coordinates().start(),
                                    variants.coordinates().representativeTranscript(),
                                    variants.coordinates().ensemblVersion(),
                                    variants.coordinates().chromosome(),
                                    variants.coordinates().referenceBuild(),
                                    idVariantGroupVariants)
                            .execute();
                }

                for (CivicVariantTypes variantTypesGroup : variants.variant_types()) {
                    context.insertInto(CIVICVARIANTSGROUPSTYPES,
                            CIVICVARIANTSGROUPSTYPES.DISPLAYNAME,
                            CIVICVARIANTSGROUPSTYPES.DESCRIPTION,
                            CIVICVARIANTSGROUPSTYPES.URL,
                            CIVICVARIANTSGROUPSTYPES.SOID,
                            CIVICVARIANTSGROUPSTYPES.IDVARIANTTYPES,
                            CIVICVARIANTSGROUPSTYPES.NAME,
                            CIVICVARIANTSGROUPSTYPES.CIVICVARIANTSGROUPSVARIANTSID)
                            .values(variantTypesGroup.displayName(),
                                    variantTypesGroup.description(),
                                    variantTypesGroup.url(),
                                    variantTypesGroup.soId(),
                                    variantTypesGroup.id(),
                                    variantTypesGroup.name(),
                                    idVariantGroupVariants)
                            .execute();
                }
            }
        }

        for (CivicEvidenceItems evidenceItems : civic.evidenceItem()) {
            int idEvidenceItems = context.insertInto(CIVICEVIDENCEITEMS,
                    CIVICEVIDENCEITEMS.STATUS,
                    CIVICEVIDENCEITEMS.RATING,
                    CIVICEVIDENCEITEMS.DRUGINTERACTIONTYPE,
                    CIVICEVIDENCEITEMS.DESCRIPTION,
                    CIVICEVIDENCEITEMS.OPENCHANGECOUNT,
                    CIVICEVIDENCEITEMS.EVIDENCETYPE,
                    CIVICEVIDENCEITEMS.VARIANTORIGIN,
                    CIVICEVIDENCEITEMS.EVIDENCEDIRECTION,
                    CIVICEVIDENCEITEMS.VARIANTID,
                    CIVICEVIDENCEITEMS.CLINICALSIGNIFICANCE,
                    CIVICEVIDENCEITEMS.EVIDENCELEVEL,
                    CIVICEVIDENCEITEMS.TYPE,
                    CIVICEVIDENCEITEMS.IDEVIDENCEITEMS,
                    CIVICEVIDENCEITEMS.NAME,
                    CIVICEVIDENCEITEMS.CIVICID)
                    .values(evidenceItems.status(),
                            evidenceItems.rating(),
                            evidenceItems.drugInteractionType(),
                            evidenceItems.description(),
                            evidenceItems.openChangeCount(),
                            evidenceItems.evidenceType(),
                            evidenceItems.variantOrigin(),
                            evidenceItems.evidenceDirection(),
                            evidenceItems.variantId(),
                            evidenceItems.clinicalSignificance(),
                            evidenceItems.evidenceLevel(),
                            evidenceItems.type(),
                            evidenceItems.id(),
                            evidenceItems.name(),
                            id)
                    .returning(CIVICEVIDENCEITEMS.ID)
                    .fetchOne()
                    .getValue(CIVICEVIDENCEITEMS.ID);

            for (CivicDrugs drugs : evidenceItems.drugs()) {
                context.insertInto(CIVICDRUGS, CIVICDRUGS.PUBCHEMID, CIVICDRUGS.IDDRUGS, CIVICDRUGS.NAME, CIVICDRUGS.CIVICEVIDENCEITEMSID)
                        .values(drugs.pubchemId(), drugs.id(), drugs.name(), idEvidenceItems)
                        .execute();
            }

            context.insertInto(CIVICDISEASE,
                    CIVICDISEASE.DOID,
                    CIVICDISEASE.URL,
                    CIVICDISEASE.DISPLAYNAME,
                    CIVICDISEASE.IDDISEASE,
                    CIVICDISEASE.NAME,
                    CIVICDISEASE.CIVICEVIDENCEITEMSID)
                    .values(evidenceItems.disease().doid(),
                            evidenceItems.disease().url(),
                            evidenceItems.disease().displayName(),
                            evidenceItems.disease().id(),
                            evidenceItems.disease().name(),
                            idEvidenceItems)
                    .execute();

            int idEvidenceItemsSource = context.insertInto(CIVICEVIDENCEITEMSSOURCE,
                    CIVICEVIDENCEITEMSSOURCE.STATUS,
                    CIVICEVIDENCEITEMSSOURCE.OPENACCESS,
                    CIVICEVIDENCEITEMSSOURCE.NAME,
                    CIVICEVIDENCEITEMSSOURCE.JOURNAL,
                    CIVICEVIDENCEITEMSSOURCE.CITATION,
                    CIVICEVIDENCEITEMSSOURCE.PMC_ID,
                    CIVICEVIDENCEITEMSSOURCE.FULLJOURNALTITLE,
                    CIVICEVIDENCEITEMSSOURCE.SOURCEURL,
                    CIVICEVIDENCEITEMSSOURCE.PUBMEDID,
                    CIVICEVIDENCEITEMSSOURCE.ISREVIEW,
                    CIVICEVIDENCEITEMSSOURCE.IDSOURCE,
                    CIVICEVIDENCEITEMSSOURCE.CIVICEVIDENCEITEMSID)
                    .values(evidenceItems.source().status(),
                            evidenceItems.source().openAccess(),
                            evidenceItems.source().name(),
                            evidenceItems.source().journal(),
                            evidenceItems.source().citation(),
                            evidenceItems.source().pmc_Id(),
                            evidenceItems.source().fullJournalTitle(),
                            evidenceItems.source().sourceUrl(),
                            evidenceItems.source().pubmedId(),
                            evidenceItems.source().isReview(),
                            evidenceItems.source().id(),
                            idEvidenceItems)
                    .returning(CIVICEVIDENCEITEMSSOURCE.ID)
                    .fetchOne()
                    .getValue(CIVICEVIDENCEITEMSSOURCE.ID);

            context.insertInto(CIVICEVIDENCEITEMSPUBLICATION,
                    CIVICEVIDENCEITEMSPUBLICATION.YEAR,
                    CIVICEVIDENCEITEMSPUBLICATION.DAY,
                    CIVICEVIDENCEITEMSPUBLICATION.MONTH,
                    CIVICEVIDENCEITEMSPUBLICATION.CIVICEVIDENCEITEMSSOURCEID)
                    .values(evidenceItems.source().publicationDate().year(),
                            evidenceItems.source().publicationDate().day(),
                            evidenceItems.source().publicationDate().month(),
                            idEvidenceItemsSource)
                    .execute();

            for (CivicClinicalTrial clinicalTrial : evidenceItems.source().clinicalTrials()) {
                context.insertInto(CIVICEVIDENCEITEMSCLINICALTRIAL,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.NCT_ID,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.DESCRIPTION,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.CLINICAL_TRIAL_URL,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.NAME,
                        CIVICEVIDENCEITEMSCLINICALTRIAL.CIVICEVIDENCEITEMSSOURCEID)
                        .values(clinicalTrial.nct_id(),
                                clinicalTrial.description(),
                                clinicalTrial.clinical_trial_url(),
                                clinicalTrial.name(),
                                idEvidenceItemsSource)
                        .execute();
            }
        }

        if (civic.sources() != null) {
            for (CivicSource source : civic.sources()) {
                int idSource = context.insertInto(CIVICSOURCE,
                        CIVICSOURCE.STATUS,
                        CIVICSOURCE.OPENACCESS,
                        CIVICSOURCE.NAME,
                        CIVICSOURCE.JOURNAL,
                        CIVICSOURCE.CITATION,
                        CIVICSOURCE.PMC_ID,
                        CIVICSOURCE.FULLJOURNALTITLE,
                        CIVICSOURCE.SOURCEURL,
                        CIVICSOURCE.PUBMEDID,
                        CIVICSOURCE.ISREVIEW,
                        CIVICSOURCE.IDSOURCE,
                        CIVICSOURCE.CIVICID)
                        .values(source.status(),
                                source.openAccess(),
                                source.name(),
                                source.journal(),
                                source.citation(),
                                source.pmc_Id(),
                                source.fullJournalTitle(),
                                source.sourceUrl(),
                                source.pubmedId(),
                                source.isReview(),
                                source.id(),
                                id)
                        .returning(CIVICSOURCE.ID)
                        .fetchOne()
                        .getValue(CIVICSOURCE.ID);

                context.insertInto(CIVICPUBLICATION,
                        CIVICPUBLICATION.YEAR,
                        CIVICPUBLICATION.DAY,
                        CIVICPUBLICATION.MONTH,
                        CIVICPUBLICATION.CIVICSOURCEID)
                        .values(source.publicationDate().year(), source.publicationDate().day(), source.publicationDate().month(), idSource)
                        .execute();

                for (CivicClinicalTrial clinicalTrial : source.clinicalTrials()) {
                    context.insertInto(CIVICCLINICALTRIAL,
                            CIVICCLINICALTRIAL.NCT_ID,
                            CIVICCLINICALTRIAL.DESCRIPTION,
                            CIVICCLINICALTRIAL.CLINICAL_TRIAL_URL,
                            CIVICCLINICALTRIAL.NAME,
                            CIVICCLINICALTRIAL.CIVICSOURCEID)
                            .values(clinicalTrial.nct_id(),
                                    clinicalTrial.description(),
                                    clinicalTrial.clinical_trial_url(),
                                    clinicalTrial.name(),
                                    idSource)
                            .execute();
                }
            }
        }

        context.insertInto(CIVICERROR, CIVICERROR.CIVICID).values(id).execute();

        int idLifeActions = context.insertInto(CIVICLIFECYCLEACTIONS, CIVICLIFECYCLEACTIONS.CIVICID)
                .values(id)
                .returning(CIVICLIFECYCLEACTIONS.ID)
                .fetchOne()
                .getValue(CIVICLIFECYCLEACTIONS.ID);

        if (civic.lifecycleActions().lastCommentedOn() != null) {
            int idLastCommendOn =
                    context.insertInto(CIVICLASTCOMMENTEDON, CIVICLASTCOMMENTEDON.TIMESTAMP, CIVICLASTCOMMENTEDON.CIVICLIFECYCLEACTIONSID)
                            .values(civic.lifecycleActions().lastCommentedOn().timestamp(), idLifeActions)
                            .returning(CIVICLASTCOMMENTEDON.ID)
                            .fetchOne()
                            .getValue(CIVICLASTCOMMENTEDON.ID);

            CivicUser userLastCommend = civic.lifecycleActions().lastCommentedOn().user();
            int idLastCommentUser = context.insertInto(CIVICLASTCOMMENTEDONUSER,
                    CIVICLASTCOMMENTEDONUSER.USERNAME,
                    CIVICLASTCOMMENTEDONUSER.AREAOFEXPERTISE,
                    CIVICLASTCOMMENTEDONUSER.TWITTERHANDLE,
                    CIVICLASTCOMMENTEDONUSER.NAME,
                    CIVICLASTCOMMENTEDONUSER.BIO,
                    CIVICLASTCOMMENTEDONUSER.URL,
                    CIVICLASTCOMMENTEDONUSER.CREATEDAT,
                    CIVICLASTCOMMENTEDONUSER.ACCEPTEDLICENSE,
                    CIVICLASTCOMMENTEDONUSER.AFFILIATION,
                    CIVICLASTCOMMENTEDONUSER.AVATARURL,
                    CIVICLASTCOMMENTEDONUSER.ROLE,
                    CIVICLASTCOMMENTEDONUSER.FACEBOOKPROFILE,
                    CIVICLASTCOMMENTEDONUSER.LINKEDINPROFILE,
                    CIVICLASTCOMMENTEDONUSER.ORCID,
                    CIVICLASTCOMMENTEDONUSER.DISPLAYNAME,
                    CIVICLASTCOMMENTEDONUSER.LASTSEENAT,
                    CIVICLASTCOMMENTEDONUSER.FEATUREDEXPERT,
                    CIVICLASTCOMMENTEDONUSER.IDUSER,
                    CIVICLASTCOMMENTEDONUSER.SIGNUPCOMPLETE,
                    CIVICLASTCOMMENTEDONUSER.CIVICLASTCOMMENTEDONID)
                    .values(userLastCommend.username(),
                            userLastCommend.areaOfExpertise(),
                            userLastCommend.twitterHandle(),
                            userLastCommend.name(),
                            userLastCommend.bio(),
                            userLastCommend.url(),
                            userLastCommend.createdAt(),
                            userLastCommend.acceptedLicense(),
                            userLastCommend.affiliation(),
                            userLastCommend.avatarUrl(),
                            userLastCommend.role(),
                            userLastCommend.facebookProfile(),
                            userLastCommend.linkedinProfile(),
                            userLastCommend.orcid(),
                            userLastCommend.displayName(),
                            userLastCommend.lastSeenAt(),
                            userLastCommend.featuredExpert(),
                            userLastCommend.id(),
                            userLastCommend.signupComplete(),
                            idLastCommendOn)
                    .returning(CIVICLASTCOMMENTEDONUSER.ID)
                    .fetchOne()
                    .getValue(CIVICLASTCOMMENTEDONUSER.ID);

            context.insertInto(CIVICLASTCOMMENTEDONAVATARS,
                    CIVICLASTCOMMENTEDONAVATARS.X32,
                    CIVICLASTCOMMENTEDONAVATARS.X14,
                    CIVICLASTCOMMENTEDONAVATARS.X64,
                    CIVICLASTCOMMENTEDONAVATARS.X128,
                    CIVICLASTCOMMENTEDONAVATARS.CIVICLASTCOMMENTEDONUSERID)
                    .values(userLastCommend.avatars().x32(),
                            userLastCommend.avatars().x14(),
                            userLastCommend.avatars().x64(),
                            userLastCommend.avatars().x128(),
                            idLastCommentUser)
                    .execute();

            int idLastCommentOnOrganization = context.insertInto(CIVICLASTCOMMENTEDONORGANIZATION,
                    CIVICLASTCOMMENTEDONORGANIZATION.URL,
                    CIVICLASTCOMMENTEDONORGANIZATION.IDORGANIZATION,
                    CIVICLASTCOMMENTEDONORGANIZATION.DESCRIPTION,
                    CIVICLASTCOMMENTEDONORGANIZATION.NAME,
                    CIVICLASTCOMMENTEDONORGANIZATION.CIVICLASTCOMMENTEDONUSERID)
                    .values(userLastCommend.organization().url(),
                            userLastCommend.organization().id(),
                            userLastCommend.organization().description(),
                            userLastCommend.organization().name(),
                            idLastCommentUser)
                    .returning(CIVICLASTCOMMENTEDONORGANIZATION.ID)
                    .fetchOne()
                    .getValue(CIVICLASTCOMMENTEDONORGANIZATION.ID);

            context.insertInto(CIVICLASTCOMMENTEDONPROFILEIMAGE,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X32,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X256,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X14,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X64,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.X128,
                    CIVICLASTCOMMENTEDONPROFILEIMAGE.CIVICLASTCOMMENTEDONORGANIZATIONID)
                    .values(userLastCommend.organization().profileImage().x32(),
                            userLastCommend.organization().profileImage().x256(),
                            userLastCommend.organization().profileImage().x14(),
                            userLastCommend.organization().profileImage().x64(),
                            userLastCommend.organization().profileImage().x128(),
                            idLastCommentOnOrganization)
                    .execute();
        }

        if (civic.lifecycleActions().lastModified() != null) {
            int idLastModiefied =
                    context.insertInto(CIVICLASTMODIFIED, CIVICLASTMODIFIED.TIMESTAMP, CIVICLASTMODIFIED.CIVICLIFECYCLEACTIONSID)
                            .values(civic.lifecycleActions().lastModified().timestamp(), idLifeActions)
                            .returning(CIVICLASTMODIFIED.ID)
                            .fetchOne()
                            .getValue(CIVICLASTMODIFIED.ID);

            CivicUser userModified = civic.lifecycleActions().lastModified().user();
            int idLastModiefiedUser = context.insertInto(CIVICLASTMODIFIEDUSER,
                    CIVICLASTMODIFIEDUSER.USERNAME,
                    CIVICLASTMODIFIEDUSER.AREAOFEXPERTISE,
                    CIVICLASTMODIFIEDUSER.TWITTERHANDLE,
                    CIVICLASTMODIFIEDUSER.NAME,
                    CIVICLASTMODIFIEDUSER.BIO,
                    CIVICLASTMODIFIEDUSER.URL,
                    CIVICLASTMODIFIEDUSER.CREATEDAT,
                    CIVICLASTMODIFIEDUSER.ACCEPTEDLICENSE,
                    CIVICLASTMODIFIEDUSER.AFFILIATION,
                    CIVICLASTMODIFIEDUSER.AVATARURL,
                    CIVICLASTMODIFIEDUSER.ROLE,
                    CIVICLASTMODIFIEDUSER.FACEBOOKPROFILE,
                    CIVICLASTMODIFIEDUSER.LINKEDINPROFILE,
                    CIVICLASTMODIFIEDUSER.ORCID,
                    CIVICLASTMODIFIEDUSER.DISPLAYNAME,
                    CIVICLASTMODIFIEDUSER.LASTSEENAT,
                    CIVICLASTMODIFIEDUSER.FEATUREDEXPERT,
                    CIVICLASTMODIFIEDUSER.IDUSER,
                    CIVICLASTMODIFIEDUSER.SIGNUPCOMPLETE,
                    CIVICLASTMODIFIEDUSER.CIVICLASTMODIFIEDID)
                    .values(userModified.username(),
                            userModified.areaOfExpertise(),
                            userModified.twitterHandle(),
                            userModified.name(),
                            userModified.bio(),
                            userModified.url(),
                            userModified.createdAt(),
                            userModified.acceptedLicense(),
                            userModified.affiliation(),
                            userModified.avatarUrl(),
                            userModified.role(),
                            userModified.facebookProfile(),
                            userModified.linkedinProfile(),
                            userModified.orcid(),
                            userModified.displayName(),
                            userModified.lastSeenAt(),
                            userModified.featuredExpert(),
                            userModified.id(),
                            userModified.signupComplete(),
                            idLastModiefied)
                    .returning(CIVICLASTMODIFIEDUSER.ID)
                    .fetchOne()
                    .getValue(CIVICLASTMODIFIEDUSER.ID);

            context.insertInto(CIVICLASTMODIFIEDAVATARS,
                    CIVICLASTMODIFIEDAVATARS.X32,
                    CIVICLASTMODIFIEDAVATARS.X14,
                    CIVICLASTMODIFIEDAVATARS.X64,
                    CIVICLASTMODIFIEDAVATARS.X128,
                    CIVICLASTMODIFIEDAVATARS.CIVICLASTMODIFIEDID)
                    .values(userModified.avatars().x32(),
                            userModified.avatars().x14(),
                            userModified.avatars().x64(),
                            userModified.avatars().x128(),
                            idLastModiefiedUser)
                    .execute();

            int idLastModiefiedOrganization = context.insertInto(CIVICLASTMODIFIEDORGANIZATION,
                    CIVICLASTMODIFIEDORGANIZATION.URL,
                    CIVICLASTMODIFIEDORGANIZATION.IDORGANIZATION,
                    CIVICLASTMODIFIEDORGANIZATION.DESCRIPTION,
                    CIVICLASTMODIFIEDORGANIZATION.NAME,
                    CIVICLASTMODIFIEDORGANIZATION.CIVICLASTMODIFIEDUSERID)
                    .values(userModified.organization().url(),
                            userModified.organization().id(),
                            userModified.organization().description(),
                            userModified.organization().name(),
                            idLastModiefiedUser)
                    .returning(CIVICLASTMODIFIEDORGANIZATION.ID)
                    .fetchOne()
                    .getValue(CIVICLASTMODIFIEDORGANIZATION.ID);

            if (userModified.organization().profileImage() != null) {
                context.insertInto(CIVICLASTMODIFIEDPROFILEIMAGE,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X32,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X256,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X14,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X64,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X128,
                        CIVICLASTMODIFIEDPROFILEIMAGE.CIVICLASTMODIFIEDORGANIZATIONID)
                        .values(userModified.organization().profileImage().x32(),
                                userModified.organization().profileImage().x256(),
                                userModified.organization().profileImage().x14(),
                                userModified.organization().profileImage().x64(),
                                userModified.organization().profileImage().x128(),
                                idLastModiefiedOrganization)
                        .execute();
            }

        }

        if (civic.lifecycleActions().lastReviewed() != null) {
            int idLastReviewed =
                    context.insertInto(CIVICLASTREVIEWED, CIVICLASTREVIEWED.TIMESTAMP, CIVICLASTREVIEWED.CIVICLIFECYCLEACTIONSID)
                            .values(civic.lifecycleActions().lastCommentedOn().timestamp(), idLifeActions)
                            .returning(CIVICLASTREVIEWED.ID)
                            .fetchOne()
                            .getValue(CIVICLASTREVIEWED.ID);

            CivicUser userLastReviewed = civic.lifecycleActions().lastReviewed().user();
            int idLastReviewedUser = context.insertInto(CIVICLASTREVIEWEDUSER,
                    CIVICLASTREVIEWEDUSER.USERNAME,
                    CIVICLASTREVIEWEDUSER.AREAOFEXPERTISE,
                    CIVICLASTREVIEWEDUSER.TWITTERHANDLE,
                    CIVICLASTREVIEWEDUSER.NAME,
                    CIVICLASTREVIEWEDUSER.BIO,
                    CIVICLASTREVIEWEDUSER.URL,
                    CIVICLASTREVIEWEDUSER.CREATEDAT,
                    CIVICLASTREVIEWEDUSER.ACCEPTEDLICENSE,
                    CIVICLASTREVIEWEDUSER.AFFILIATION,
                    CIVICLASTREVIEWEDUSER.AVATARURL,
                    CIVICLASTREVIEWEDUSER.ROLE,
                    CIVICLASTREVIEWEDUSER.FACEBOOKPROFILE,
                    CIVICLASTREVIEWEDUSER.LINKEDINPROFILE,
                    CIVICLASTREVIEWEDUSER.ORCID,
                    CIVICLASTREVIEWEDUSER.DISPLAYNAME,
                    CIVICLASTREVIEWEDUSER.LASTSEENAT,
                    CIVICLASTREVIEWEDUSER.FEATUREDEXPERT,
                    CIVICLASTREVIEWEDUSER.IDUSER,
                    CIVICLASTREVIEWEDUSER.SIGNUPCOMPLETE,
                    CIVICLASTREVIEWEDUSER.CIVICLASTREVIEWEDID)
                    .values(userLastReviewed.username(),
                            userLastReviewed.areaOfExpertise(),
                            userLastReviewed.twitterHandle(),
                            userLastReviewed.name(),
                            userLastReviewed.bio(),
                            userLastReviewed.url(),
                            userLastReviewed.createdAt(),
                            userLastReviewed.acceptedLicense(),
                            userLastReviewed.affiliation(),
                            userLastReviewed.avatarUrl(),
                            userLastReviewed.role(),
                            userLastReviewed.facebookProfile(),
                            userLastReviewed.linkedinProfile(),
                            userLastReviewed.orcid(),
                            userLastReviewed.displayName(),
                            userLastReviewed.lastSeenAt(),
                            userLastReviewed.featuredExpert(),
                            userLastReviewed.id(),
                            userLastReviewed.signupComplete(),
                            idLastReviewed)
                    .returning(CIVICLASTREVIEWEDUSER.ID)
                    .fetchOne()
                    .getValue(CIVICLASTREVIEWEDUSER.ID);

            context.insertInto(CIVICLASTREVIEWEDAVATARS,
                    CIVICLASTREVIEWEDAVATARS.X32,
                    CIVICLASTREVIEWEDAVATARS.X14,
                    CIVICLASTREVIEWEDAVATARS.X64,
                    CIVICLASTREVIEWEDAVATARS.X128,
                    CIVICLASTREVIEWEDAVATARS.CIVICLASTREVIEWEDID)
                    .values(userLastReviewed.avatars().x32(),
                            userLastReviewed.avatars().x14(),
                            userLastReviewed.avatars().x64(),
                            userLastReviewed.avatars().x128(),
                            idLastReviewedUser)
                    .execute();

            int idLastReviewedOrganization = context.insertInto(CIVICLASTREVIEWEDORGANIZATION,
                    CIVICLASTREVIEWEDORGANIZATION.URL,
                    CIVICLASTREVIEWEDORGANIZATION.IDORGANIZATION,
                    CIVICLASTREVIEWEDORGANIZATION.DESCRIPTION,
                    CIVICLASTREVIEWEDORGANIZATION.NAME,
                    CIVICLASTREVIEWEDORGANIZATION.CIVICLASTREVIEWEDUSERID)
                    .values(userLastReviewed.organization().url(),
                            userLastReviewed.organization().id(),
                            userLastReviewed.organization().description(),
                            userLastReviewed.organization().name(),
                            idLastReviewedUser)
                    .returning(CIVICLASTREVIEWEDORGANIZATION.ID)
                    .fetchOne()
                    .getValue(CIVICLASTREVIEWEDORGANIZATION.ID);

            context.insertInto(CIVICLASTREVIEWEDPROFILEIMAGE,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X32,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X256,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X14,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X64,
                    CIVICLASTREVIEWEDPROFILEIMAGE.X128,
                    CIVICLASTREVIEWEDPROFILEIMAGE.CIVICLASTREVIEWEDORGANIZATIONID)
                    .values(userLastReviewed.organization().profileImage().x32(),
                            userLastReviewed.organization().profileImage().x256(),
                            userLastReviewed.organization().profileImage().x14(),
                            userLastReviewed.organization().profileImage().x64(),
                            userLastReviewed.organization().profileImage().x128(),
                            idLastReviewedOrganization)
                    .execute();
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        context.deleteFrom(CIVIC).execute();
        context.deleteFrom(CIVICASSERTIONS).execute();
        context.deleteFrom(CIVICHGVSEXPRESSIONS).execute();
        context.deleteFrom(CIVICCLINVARENTRIES).execute();
        context.deleteFrom(CIVICVARIANTALIASES).execute();
        context.deleteFrom(CIVICVARIANTTYPES).execute();
        context.deleteFrom(CIVICDESCRIPTION).execute();
        context.deleteFrom(CIVICCOORDINATES).execute();
        context.deleteFrom(CIVICVARIANTSGROUPS).execute();
        context.deleteFrom(CIVICVARIANTSGROUPSVARIANTS).execute();
        context.deleteFrom(CIVICVARIANTSGROUPSCOORDINATES).execute();
        context.deleteFrom(CIVICVARIANTSGROUPSTYPES).execute();
        context.deleteFrom(CIVICEVIDENCEITEMS).execute();
        context.deleteFrom(CIVICDRUGS).execute();
        context.deleteFrom(CIVICDISEASE).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSSOURCE).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSPUBLICATION).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSCLINICALTRIAL).execute();
        context.deleteFrom(CIVICSOURCE).execute();
        context.deleteFrom(CIVICERROR).execute();
        context.deleteFrom(CIVICPUBLICATION).execute();
        context.deleteFrom(CIVICCLINICALTRIAL).execute();
        context.deleteFrom(CIVICLIFECYCLEACTIONS).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDON).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONUSER).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONAVATARS).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONORGANIZATION).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONPROFILEIMAGE).execute();
        context.deleteFrom(CIVICLASTMODIFIED).execute();
        context.deleteFrom(CIVICLASTMODIFIEDUSER).execute();
        context.deleteFrom(CIVICLASTMODIFIEDAVATARS).execute();
        context.deleteFrom(CIVICLASTMODIFIEDORGANIZATION).execute();
        context.deleteFrom(CIVICLASTMODIFIEDPROFILEIMAGE).execute();
        context.deleteFrom(CIVICLASTREVIEWED).execute();
        context.deleteFrom(CIVICLASTREVIEWEDUSER).execute();
        context.deleteFrom(CIVICLASTREVIEWEDAVATARS).execute();
        context.deleteFrom(CIVICLASTREVIEWEDORGANIZATION).execute();
        context.deleteFrom(CIVICLASTREVIEWEDPROFILEIMAGE).execute();
    }
}

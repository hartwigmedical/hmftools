package com.hartwig.hmftools.vicc.dao;

import static com.hartwig.hmftools.vicc.database.Tables.CIVIC;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICASSERTION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCLINICALTRIAL;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCLINVARENTRY;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDISEASE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICDRUG;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEM;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMCLINICALTRIAL;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMPUBLICATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICEVIDENCEITEMSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICHGVSEXPRESSION;
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
import static com.hartwig.hmftools.vicc.database.Tables.CIVICPROVISIONALVALUE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICPUBLICATION;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICSOURCE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTALIAS;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTGROUP;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTGROUPCOORDINATES;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTGROUPTYPE;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTGROUPVARIANT;
import static com.hartwig.hmftools.vicc.database.Tables.CIVICVARIANTTYPE;

import com.hartwig.hmftools.vicc.datamodel.civic.Civic;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicClinicalTrial;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicCoordinates;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicDrug;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicLastCommentedOn;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicLastModified;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicLastReviewed;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicProfileImage;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicProvisionalValue;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicSource;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicUser;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicVariant;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicVariantGroup;
import com.hartwig.hmftools.vicc.datamodel.civic.CivicVariantType;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

final class CivicDAOFunctions {

    private CivicDAOFunctions() {
    }

    static void write(@NotNull DSLContext context, int viccEntryId, @NotNull Civic civic) {
        int id = context.insertInto(CIVIC,
                CIVIC.ENTREZID,
                CIVIC.ENTREZNAME,
                CIVIC.NAME,
                CIVIC.TYPE,
                CIVIC.CIVICACTIONABILITYSCORE,
                CIVIC.ALLELEREGISTRYID,
                CIVIC.IDCIVIC,
                CIVIC.GENEID,
                CIVIC.DESCRIPTION,
                CIVIC.VICCENTRYID)
                .values(civic.entrezId(),
                        civic.entrezName(),
                        civic.name(),
                        civic.type(),
                        civic.civicActionabilityScore(),
                        civic.alleleRegistryId(),
                        civic.id(),
                        civic.geneId(),
                        civic.description(),
                        viccEntryId)
                .returning(CIVIC.ID)
                .fetchOne()
                .getValue(CIVIC.ID);

        for (String assertion : civic.assertions()) {
            context.insertInto(CIVICASSERTION, CIVICASSERTION.ASSERTION, CIVICASSERTION.CIVICID).values(assertion, id).execute();
        }

        for (String hgvsExpression : civic.hgvsExpressions()) {
            context.insertInto(CIVICHGVSEXPRESSION, CIVICHGVSEXPRESSION.HGVSEXPRESSION, CIVICHGVSEXPRESSION.CIVICID)
                    .values(hgvsExpression, id)
                    .execute();
        }

        for (String clinVarEntry : civic.clinVarEntries()) {
            context.insertInto(CIVICCLINVARENTRY, CIVICCLINVARENTRY.CLINVARENTRY, CIVICCLINVARENTRY.CIVICID)
                    .values(clinVarEntry, id)
                    .execute();
        }

        for (String variantAlias : civic.variantAliases()) {
            context.insertInto(CIVICVARIANTALIAS, CIVICVARIANTALIAS.VARIANTALIAS, CIVICVARIANTALIAS.CIVICID)
                    .values(variantAlias, id)
                    .execute();
        }

        for (CivicVariantType variantType : civic.variantTypes()) {
            context.insertInto(CIVICVARIANTTYPE,
                    CIVICVARIANTTYPE.NAME,
                    CIVICVARIANTTYPE.DISPLAYNAME,
                    CIVICVARIANTTYPE.DESCRIPTION,
                    CIVICVARIANTTYPE.URL,
                    CIVICVARIANTTYPE.SOID,
                    CIVICVARIANTTYPE.IDVARIANTTYPE,
                    CIVICVARIANTTYPE.CIVICID)
                    .values(variantType.name(),
                            variantType.displayName(),
                            variantType.description(),
                            variantType.url(),
                            variantType.soId(),
                            variantType.id(),
                            id)
                    .execute();
        }

        CivicProvisionalValue provisionalValue = civic.provisionalValue();
        if (provisionalValue != null) {
            context.insertInto(CIVICPROVISIONALVALUE,
                    CIVICPROVISIONALVALUE.REVISIONID,
                    CIVICPROVISIONALVALUE.VALUE,
                    CIVICPROVISIONALVALUE.CIVICID).values(provisionalValue.revisionId(), provisionalValue.value(), id).execute();
        }

        context.insertInto(CIVICCOORDINATES,
                CIVICCOORDINATES.CHROMOSOME,
                CIVICCOORDINATES.START,
                CIVICCOORDINATES.STOP,
                CIVICCOORDINATES.REFERENCEBASES,
                CIVICCOORDINATES.VARIANTBASES,
                CIVICCOORDINATES.REPRESENTATIVETRANSCRIPT,
                CIVICCOORDINATES.ENSEMBLVERSION,
                CIVICCOORDINATES.REFERENCEBUILD,
                CIVICCOORDINATES.CHROMOSOME2,
                CIVICCOORDINATES.START2,
                CIVICCOORDINATES.STOP2,
                CIVICCOORDINATES.REPRESENTATIVETRANSCRIPT2,
                CIVICCOORDINATES.CIVICID)
                .values(civic.coordinates().chromosome(),
                        civic.coordinates().start(),
                        civic.coordinates().stop(),
                        civic.coordinates().referenceBases(),
                        civic.coordinates().variantBases(),
                        civic.coordinates().representativeTranscript(),
                        civic.coordinates().ensemblVersion(),
                        civic.coordinates().referenceBuild(),
                        civic.coordinates().chromosome2(),
                        civic.coordinates().start2(),
                        civic.coordinates().stop2(),
                        civic.coordinates().representativeTranscript2(),
                        id)
                .execute();

        for (CivicVariantGroup variantGroup : civic.variantGroups()) {
            int idVariantGroup = context.insertInto(CIVICVARIANTGROUP,
                    CIVICVARIANTGROUP.NAME,
                    CIVICVARIANTGROUP.TYPE,
                    CIVICVARIANTGROUP.DESCRIPTION,
                    CIVICVARIANTGROUP.IDVARIANTGROUP,
                    CIVICVARIANTGROUP.CIVICID)
                    .values(variantGroup.name(), variantGroup.type(), variantGroup.description(), variantGroup.id(), id)
                    .returning(CIVICVARIANTGROUP.ID)
                    .fetchOne()
                    .getValue(CIVICVARIANTGROUP.ID);

            for (CivicVariant variant : variantGroup.variants()) {
                int idVariantGroupVariant = context.insertInto(CIVICVARIANTGROUPVARIANT,
                        CIVICVARIANTGROUPVARIANT.ENTREZID,
                        CIVICVARIANTGROUPVARIANT.ENTREZNAME,
                        CIVICVARIANTGROUPVARIANT.NAME,
                        CIVICVARIANTGROUPVARIANT.TYPE,
                        CIVICVARIANTGROUPVARIANT.CIVICACTIONABILITYSCORE,
                        CIVICVARIANTGROUPVARIANT.IDVARIANT,
                        CIVICVARIANTGROUPVARIANT.GENEID,
                        CIVICVARIANTGROUPVARIANT.DESCRIPTION,
                        CIVICVARIANTGROUPVARIANT.CIVICVARIANTGROUPID)
                        .values(variant.entrezId(),
                                variant.entrezName(),
                                variant.name(),
                                variant.type(),
                                variant.civicActionabilityScore(),
                                variant.id(),
                                variant.geneId(),
                                variant.description(),
                                idVariantGroup)
                        .returning(CIVICVARIANTGROUPVARIANT.ID)
                        .fetchOne()
                        .getValue(CIVICVARIANTGROUPVARIANT.ID);

                CivicCoordinates coordinates = variant.coordinates();
                if (coordinates != null) {
                    context.insertInto(CIVICVARIANTGROUPCOORDINATES,
                            CIVICVARIANTGROUPCOORDINATES.CHROMOSOME,
                            CIVICVARIANTGROUPCOORDINATES.START,
                            CIVICVARIANTGROUPCOORDINATES.STOP,
                            CIVICVARIANTGROUPCOORDINATES.REFERENCEBASES,
                            CIVICVARIANTGROUPCOORDINATES.VARIANTBASES,
                            CIVICVARIANTGROUPCOORDINATES.REPRESENTATIVETRANSCRIPT,
                            CIVICVARIANTGROUPCOORDINATES.ENSEMBLVERSION,
                            CIVICVARIANTGROUPCOORDINATES.REFERENCEBUILD,
                            CIVICVARIANTGROUPCOORDINATES.CHROMOSOME2,
                            CIVICVARIANTGROUPCOORDINATES.START2,
                            CIVICVARIANTGROUPCOORDINATES.STOP2,
                            CIVICVARIANTGROUPCOORDINATES.REPRESENTATIVETRANSCRIPT2,
                            CIVICVARIANTGROUPCOORDINATES.CIVICVARIANTGROUPVARIANTID)
                            .values(coordinates.chromosome(),
                                    coordinates.start(),
                                    coordinates.stop(),
                                    coordinates.referenceBases(),
                                    coordinates.variantBases(),
                                    coordinates.representativeTranscript(),
                                    coordinates.ensemblVersion(),
                                    coordinates.referenceBuild(),
                                    coordinates.chromosome2(),
                                    coordinates.start2(),
                                    coordinates.stop2(),
                                    coordinates.representativeTranscript2(),
                                    idVariantGroupVariant)
                            .execute();
                }

                for (CivicVariantType variantType : variant.variantTypes()) {
                    context.insertInto(CIVICVARIANTGROUPTYPE,
                            CIVICVARIANTGROUPTYPE.NAME,
                            CIVICVARIANTGROUPTYPE.DISPLAYNAME,
                            CIVICVARIANTGROUPTYPE.DESCRIPTION,
                            CIVICVARIANTGROUPTYPE.URL,
                            CIVICVARIANTGROUPTYPE.SOID,
                            CIVICVARIANTGROUPTYPE.IDVARIANTTYPE,
                            CIVICVARIANTGROUPTYPE.CIVICVARIANTGROUPVARIANTID)
                            .values(variantType.name(),
                                    variantType.displayName(),
                                    variantType.description(),
                                    variantType.url(),
                                    variantType.soId(),
                                    variantType.id(),
                                    idVariantGroupVariant)
                            .execute();
                }
            }
        }

        int idEvidenceItem = context.insertInto(CIVICEVIDENCEITEM,
                CIVICEVIDENCEITEM.NAME,
                CIVICEVIDENCEITEM.TYPE,
                CIVICEVIDENCEITEM.STATUS,
                CIVICEVIDENCEITEM.RATING,
                CIVICEVIDENCEITEM.EVIDENCETYPE,
                CIVICEVIDENCEITEM.EVIDENCELEVEL,
                CIVICEVIDENCEITEM.EVIDENCEDIRECTION,
                CIVICEVIDENCEITEM.DRUGINTERACTIONTYPE,
                CIVICEVIDENCEITEM.VARIANTORIGIN,
                CIVICEVIDENCEITEM.CLINICALSIGNIFICANCE,
                CIVICEVIDENCEITEM.OPENCHANGECOUNT,
                CIVICEVIDENCEITEM.DESCRIPTION,
                CIVICEVIDENCEITEM.VARIANTID,
                CIVICEVIDENCEITEM.IDEVIDENCEITEM,
                CIVICEVIDENCEITEM.CIVICID)
                .values(civic.evidenceItem().name(),
                        civic.evidenceItem().type(),
                        civic.evidenceItem().status(),
                        civic.evidenceItem().rating(),
                        civic.evidenceItem().evidenceType(),
                        civic.evidenceItem().evidenceLevel(),
                        civic.evidenceItem().evidenceDirection(),
                        civic.evidenceItem().drugInteractionType(),
                        civic.evidenceItem().variantOrigin(),
                        civic.evidenceItem().clinicalSignificance(),
                        civic.evidenceItem().openChangeCount(),
                        civic.evidenceItem().description(),
                        civic.evidenceItem().variantId(),
                        civic.evidenceItem().id(),
                        id)
                .returning(CIVICEVIDENCEITEM.ID)
                .fetchOne()
                .getValue(CIVICEVIDENCEITEM.ID);

        for (CivicDrug drug : civic.evidenceItem().drugs()) {
            context.insertInto(CIVICDRUG, CIVICDRUG.NAME, CIVICDRUG.PUBCHEMID, CIVICDRUG.IDDRUG, CIVICDRUG.CIVICEVIDENCEITEMID)
                    .values(drug.name(), drug.pubchemId(), drug.id(), idEvidenceItem)
                    .execute();
        }

        context.insertInto(CIVICDISEASE,
                CIVICDISEASE.NAME,
                CIVICDISEASE.DISPLAYNAME,
                CIVICDISEASE.DOID,
                CIVICDISEASE.URL,
                CIVICDISEASE.IDDISEASE,
                CIVICDISEASE.CIVICEVIDENCEITEMID)
                .values(civic.evidenceItem().disease().name(),
                        civic.evidenceItem().disease().displayName(),
                        civic.evidenceItem().disease().doid(),
                        civic.evidenceItem().disease().url(),
                        civic.evidenceItem().disease().id(),
                        idEvidenceItem)
                .execute();

        int idEvidenceItemSource = context.insertInto(CIVICEVIDENCEITEMSOURCE,
                CIVICEVIDENCEITEMSOURCE.NAME,
                CIVICEVIDENCEITEMSOURCE.STATUS,
                CIVICEVIDENCEITEMSOURCE.OPENACCESS,
                CIVICEVIDENCEITEMSOURCE.JOURNAL,
                CIVICEVIDENCEITEMSOURCE.FULLJOURNALTITLE,
                CIVICEVIDENCEITEMSOURCE.CITATION,
                CIVICEVIDENCEITEMSOURCE.PMCID,
                CIVICEVIDENCEITEMSOURCE.SOURCEURL,
                CIVICEVIDENCEITEMSOURCE.PUBMEDID,
                CIVICEVIDENCEITEMSOURCE.ISREVIEW,
                CIVICEVIDENCEITEMSOURCE.IDSOURCE,
                CIVICEVIDENCEITEMSOURCE.CIVICEVIDENCEITEMID)
                .values(civic.evidenceItem().source().name(),
                        civic.evidenceItem().source().status(),
                        civic.evidenceItem().source().openAccess(),
                        civic.evidenceItem().source().journal(),
                        civic.evidenceItem().source().fullJournalTitle(),
                        civic.evidenceItem().source().citation(),
                        civic.evidenceItem().source().pmcId(),
                        civic.evidenceItem().source().sourceUrl(),
                        civic.evidenceItem().source().pubmedId(),
                        civic.evidenceItem().source().isReview(),
                        civic.evidenceItem().source().id(),
                        idEvidenceItem)
                .returning(CIVICEVIDENCEITEMSOURCE.ID)
                .fetchOne()
                .getValue(CIVICEVIDENCEITEMSOURCE.ID);

        context.insertInto(CIVICEVIDENCEITEMPUBLICATION,
                CIVICEVIDENCEITEMPUBLICATION.YEAR,
                CIVICEVIDENCEITEMPUBLICATION.MONTH,
                CIVICEVIDENCEITEMPUBLICATION.DAY,
                CIVICEVIDENCEITEMPUBLICATION.CIVICEVIDENCEITEMSOURCEID)
                .values(civic.evidenceItem().source().publicationDate().year(),
                        civic.evidenceItem().source().publicationDate().month(),
                        civic.evidenceItem().source().publicationDate().day(),
                        idEvidenceItemSource)
                .execute();

        for (CivicClinicalTrial clinicalTrial : civic.evidenceItem().source().clinicalTrials()) {
            context.insertInto(CIVICEVIDENCEITEMCLINICALTRIAL,
                    CIVICEVIDENCEITEMCLINICALTRIAL.NAME,
                    CIVICEVIDENCEITEMCLINICALTRIAL.NCTID,
                    CIVICEVIDENCEITEMCLINICALTRIAL.CLINICALTRIALURL,
                    CIVICEVIDENCEITEMCLINICALTRIAL.DESCRIPTION,
                    CIVICEVIDENCEITEMCLINICALTRIAL.CIVICEVIDENCEITEMSOURCEID)
                    .values(clinicalTrial.name(),
                            clinicalTrial.nctId(),
                            clinicalTrial.clinicalTrialUrl(),
                            clinicalTrial.description(),
                            idEvidenceItemSource)
                    .execute();
        }

        for (CivicSource source : civic.sources()) {
            int idSource = context.insertInto(CIVICSOURCE,
                    CIVICSOURCE.NAME,
                    CIVICSOURCE.STATUS,
                    CIVICSOURCE.OPENACCESS,
                    CIVICSOURCE.JOURNAL,
                    CIVICSOURCE.FULLJOURNALTITLE,
                    CIVICSOURCE.CITATION,
                    CIVICSOURCE.PMCID,
                    CIVICSOURCE.SOURCEURL,
                    CIVICSOURCE.PUBMEDID,
                    CIVICSOURCE.ISREVIEW,
                    CIVICSOURCE.IDSOURCE,
                    CIVICSOURCE.CIVICID)
                    .values(source.name(),
                            source.status(),
                            source.openAccess(),
                            source.journal(),
                            source.fullJournalTitle(),
                            source.citation(),
                            source.pmcId(),
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
                    CIVICPUBLICATION.MONTH,
                    CIVICPUBLICATION.DAY,
                    CIVICPUBLICATION.CIVICSOURCEID)
                    .values(source.publicationDate().year(), source.publicationDate().month(), source.publicationDate().day(), idSource)
                    .execute();

            for (CivicClinicalTrial clinicalTrial : source.clinicalTrials()) {
                context.insertInto(CIVICCLINICALTRIAL,
                        CIVICCLINICALTRIAL.NAME,
                        CIVICCLINICALTRIAL.NCTID,
                        CIVICCLINICALTRIAL.CLINICALTRIALURL,
                        CIVICCLINICALTRIAL.DESCRIPTION,
                        CIVICCLINICALTRIAL.CIVICSOURCEID)
                        .values(clinicalTrial.name(),
                                clinicalTrial.nctId(),
                                clinicalTrial.clinicalTrialUrl(),
                                clinicalTrial.description(),
                                idSource)
                        .execute();
            }
        }

        int idLifecycleActions = context.insertInto(CIVICLIFECYCLEACTIONS, CIVICLIFECYCLEACTIONS.CIVICID)
                .values(id)
                .returning(CIVICLIFECYCLEACTIONS.ID)
                .fetchOne()
                .getValue(CIVICLIFECYCLEACTIONS.ID);

        CivicLastCommentedOn lastCommentedOn = civic.lifecycleActions().lastCommentedOn();
        if (lastCommentedOn != null) {
            int idLastCommentedOn =
                    context.insertInto(CIVICLASTCOMMENTEDON, CIVICLASTCOMMENTEDON.TIMESTAMP, CIVICLASTCOMMENTEDON.CIVICLIFECYCLEACTIONSID)
                            .values(lastCommentedOn.timestamp(), idLifecycleActions)
                            .returning(CIVICLASTCOMMENTEDON.ID)
                            .fetchOne()
                            .getValue(CIVICLASTCOMMENTEDON.ID);

            CivicUser userLastCommentedOn = lastCommentedOn.user();
            int idLastCommentedOnUser = context.insertInto(CIVICLASTCOMMENTEDONUSER,
                    CIVICLASTCOMMENTEDONUSER.USERNAME,
                    CIVICLASTCOMMENTEDONUSER.NAME,
                    CIVICLASTCOMMENTEDONUSER.DISPLAYNAME,
                    CIVICLASTCOMMENTEDONUSER.ROLE,
                    CIVICLASTCOMMENTEDONUSER.AFFILIATION,
                    CIVICLASTCOMMENTEDONUSER.FEATUREDEXPERT,
                    CIVICLASTCOMMENTEDONUSER.AREAOFEXPERTISE,
                    CIVICLASTCOMMENTEDONUSER.BIO,
                    CIVICLASTCOMMENTEDONUSER.URL,
                    CIVICLASTCOMMENTEDONUSER.CREATEDAT,
                    CIVICLASTCOMMENTEDONUSER.LASTSEENAT,
                    CIVICLASTCOMMENTEDONUSER.AVATARURL,
                    CIVICLASTCOMMENTEDONUSER.TWITTERHANDLE,
                    CIVICLASTCOMMENTEDONUSER.FACEBOOKPROFILE,
                    CIVICLASTCOMMENTEDONUSER.LINKEDINPROFILE,
                    CIVICLASTCOMMENTEDONUSER.ORCID,
                    CIVICLASTCOMMENTEDONUSER.SIGNUPCOMPLETE,
                    CIVICLASTCOMMENTEDONUSER.ACCEPTEDLICENSE,
                    CIVICLASTCOMMENTEDONUSER.IDUSER,
                    CIVICLASTCOMMENTEDONUSER.CIVICLASTCOMMENTEDONID)
                    .values(userLastCommentedOn.username(),
                            userLastCommentedOn.name(),
                            userLastCommentedOn.displayName(),
                            userLastCommentedOn.role(),
                            userLastCommentedOn.affiliation(),
                            userLastCommentedOn.featuredExpert(),
                            userLastCommentedOn.areaOfExpertise(),
                            userLastCommentedOn.bio(),
                            userLastCommentedOn.url(),
                            userLastCommentedOn.createdAt(),
                            userLastCommentedOn.lastSeenAt(),
                            userLastCommentedOn.avatarUrl(),
                            userLastCommentedOn.twitterHandle(),
                            userLastCommentedOn.facebookProfile(),
                            userLastCommentedOn.linkedinProfile(),
                            userLastCommentedOn.orcid(),
                            userLastCommentedOn.signupComplete(),
                            userLastCommentedOn.acceptedLicense(),
                            userLastCommentedOn.id(),
                            idLastCommentedOn)
                    .returning(CIVICLASTCOMMENTEDONUSER.ID)
                    .fetchOne()
                    .getValue(CIVICLASTCOMMENTEDONUSER.ID);

            context.insertInto(CIVICLASTCOMMENTEDONAVATARS,
                    CIVICLASTCOMMENTEDONAVATARS.X14,
                    CIVICLASTCOMMENTEDONAVATARS.X32,
                    CIVICLASTCOMMENTEDONAVATARS.X64,
                    CIVICLASTCOMMENTEDONAVATARS.X128,
                    CIVICLASTCOMMENTEDONAVATARS.CIVICLASTCOMMENTEDONUSERID)
                    .values(userLastCommentedOn.avatars().x14(),
                            userLastCommentedOn.avatars().x32(),
                            userLastCommentedOn.avatars().x64(),
                            userLastCommentedOn.avatars().x128(),
                            idLastCommentedOnUser)
                    .execute();

            int idLastCommentOnOrganization = context.insertInto(CIVICLASTCOMMENTEDONORGANIZATION,
                    CIVICLASTCOMMENTEDONORGANIZATION.NAME,
                    CIVICLASTCOMMENTEDONORGANIZATION.URL,
                    CIVICLASTCOMMENTEDONORGANIZATION.IDORGANIZATION,
                    CIVICLASTCOMMENTEDONORGANIZATION.DESCRIPTION,
                    CIVICLASTCOMMENTEDONORGANIZATION.CIVICLASTCOMMENTEDONUSERID)
                    .values(userLastCommentedOn.organization().name(),
                            userLastCommentedOn.organization().url(),
                            userLastCommentedOn.organization().id(),
                            userLastCommentedOn.organization().description(),
                            idLastCommentedOnUser)
                    .returning(CIVICLASTCOMMENTEDONORGANIZATION.ID)
                    .fetchOne()
                    .getValue(CIVICLASTCOMMENTEDONORGANIZATION.ID);

            CivicProfileImage userLastCommentedOnProfileImage = userLastCommentedOn.organization().profileImage();
            if (userLastCommentedOnProfileImage != null) {
                context.insertInto(CIVICLASTCOMMENTEDONPROFILEIMAGE,
                        CIVICLASTCOMMENTEDONPROFILEIMAGE.X14,
                        CIVICLASTCOMMENTEDONPROFILEIMAGE.X32,
                        CIVICLASTCOMMENTEDONPROFILEIMAGE.X64,
                        CIVICLASTCOMMENTEDONPROFILEIMAGE.X128,
                        CIVICLASTCOMMENTEDONPROFILEIMAGE.X256,
                        CIVICLASTCOMMENTEDONPROFILEIMAGE.CIVICLASTCOMMENTEDONORGANIZATIONID)
                        .values(userLastCommentedOnProfileImage.x14(),
                                userLastCommentedOnProfileImage.x32(),
                                userLastCommentedOnProfileImage.x64(),
                                userLastCommentedOnProfileImage.x128(),
                                userLastCommentedOnProfileImage.x256(),
                                idLastCommentOnOrganization)
                        .execute();
            }
        }

        CivicLastModified lastModified = civic.lifecycleActions().lastModified();
        if (lastModified != null) {
            int idLastModified =
                    context.insertInto(CIVICLASTMODIFIED, CIVICLASTMODIFIED.TIMESTAMP, CIVICLASTMODIFIED.CIVICLIFECYCLEACTIONSID)
                            .values(lastModified.timestamp(), idLifecycleActions)
                            .returning(CIVICLASTMODIFIED.ID)
                            .fetchOne()
                            .getValue(CIVICLASTMODIFIED.ID);

            CivicUser userLastModified = lastModified.user();

            int idLastModifiedUser = context.insertInto(CIVICLASTMODIFIEDUSER,
                    CIVICLASTMODIFIEDUSER.USERNAME,
                    CIVICLASTMODIFIEDUSER.NAME,
                    CIVICLASTMODIFIEDUSER.DISPLAYNAME,
                    CIVICLASTMODIFIEDUSER.ROLE,
                    CIVICLASTMODIFIEDUSER.AFFILIATION,
                    CIVICLASTMODIFIEDUSER.FEATUREDEXPERT,
                    CIVICLASTMODIFIEDUSER.AREAOFEXPERTISE,
                    CIVICLASTMODIFIEDUSER.BIO,
                    CIVICLASTMODIFIEDUSER.URL,
                    CIVICLASTMODIFIEDUSER.CREATEDAT,
                    CIVICLASTMODIFIEDUSER.LASTSEENAT,
                    CIVICLASTMODIFIEDUSER.AVATARURL,
                    CIVICLASTMODIFIEDUSER.TWITTERHANDLE,
                    CIVICLASTMODIFIEDUSER.FACEBOOKPROFILE,
                    CIVICLASTMODIFIEDUSER.LINKEDINPROFILE,
                    CIVICLASTMODIFIEDUSER.ORCID,
                    CIVICLASTMODIFIEDUSER.SIGNUPCOMPLETE,
                    CIVICLASTMODIFIEDUSER.ACCEPTEDLICENSE,
                    CIVICLASTMODIFIEDUSER.IDUSER,
                    CIVICLASTMODIFIEDUSER.CIVICLASTMODIFIEDID)
                    .values(userLastModified.username(),
                            userLastModified.name(),
                            userLastModified.displayName(),
                            userLastModified.role(),
                            userLastModified.affiliation(),
                            userLastModified.featuredExpert(),
                            userLastModified.areaOfExpertise(),
                            userLastModified.bio(),
                            userLastModified.url(),
                            userLastModified.createdAt(),
                            userLastModified.lastSeenAt(),
                            userLastModified.avatarUrl(),
                            userLastModified.twitterHandle(),
                            userLastModified.facebookProfile(),
                            userLastModified.linkedinProfile(),
                            userLastModified.orcid(),
                            userLastModified.signupComplete(),
                            userLastModified.acceptedLicense(),
                            userLastModified.id(),
                            idLastModified)
                    .returning(CIVICLASTMODIFIEDUSER.ID)
                    .fetchOne()
                    .getValue(CIVICLASTMODIFIEDUSER.ID);

            context.insertInto(CIVICLASTMODIFIEDAVATARS,
                    CIVICLASTMODIFIEDAVATARS.X14,
                    CIVICLASTMODIFIEDAVATARS.X32,
                    CIVICLASTMODIFIEDAVATARS.X64,
                    CIVICLASTMODIFIEDAVATARS.X128,
                    CIVICLASTMODIFIEDAVATARS.CIVICLASTMODIFIEDID)
                    .values(userLastModified.avatars().x14(),
                            userLastModified.avatars().x32(),
                            userLastModified.avatars().x64(),
                            userLastModified.avatars().x128(),
                            idLastModifiedUser)
                    .execute();

            int idLastModifiedOrganization = context.insertInto(CIVICLASTMODIFIEDORGANIZATION,
                    CIVICLASTMODIFIEDORGANIZATION.NAME,
                    CIVICLASTMODIFIEDORGANIZATION.URL,
                    CIVICLASTMODIFIEDORGANIZATION.IDORGANIZATION,
                    CIVICLASTMODIFIEDORGANIZATION.DESCRIPTION,
                    CIVICLASTMODIFIEDORGANIZATION.CIVICLASTMODIFIEDUSERID)
                    .values(userLastModified.organization().name(),
                            userLastModified.organization().url(),
                            userLastModified.organization().id(),
                            userLastModified.organization().description(),
                            idLastModifiedUser)
                    .returning(CIVICLASTMODIFIEDORGANIZATION.ID)
                    .fetchOne()
                    .getValue(CIVICLASTMODIFIEDORGANIZATION.ID);

            CivicProfileImage userLastModifiedProfileImage = userLastModified.organization().profileImage();
            if (userLastModifiedProfileImage != null) {
                context.insertInto(CIVICLASTMODIFIEDPROFILEIMAGE,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X14,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X32,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X64,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X128,
                        CIVICLASTMODIFIEDPROFILEIMAGE.X256,
                        CIVICLASTMODIFIEDPROFILEIMAGE.CIVICLASTMODIFIEDORGANIZATIONID)
                        .values(userLastModifiedProfileImage.x14(),
                                userLastModifiedProfileImage.x32(),
                                userLastModifiedProfileImage.x64(),
                                userLastModifiedProfileImage.x128(),
                                userLastModifiedProfileImage.x256(),
                                idLastModifiedOrganization)
                        .execute();
            }
        }

        CivicLastReviewed lastReviewed = civic.lifecycleActions().lastReviewed();
        if (lastReviewed != null) {
            int idLastReviewed =
                    context.insertInto(CIVICLASTREVIEWED, CIVICLASTREVIEWED.TIMESTAMP, CIVICLASTREVIEWED.CIVICLIFECYCLEACTIONSID)
                            .values(lastReviewed.timestamp(), idLifecycleActions)
                            .returning(CIVICLASTREVIEWED.ID)
                            .fetchOne()
                            .getValue(CIVICLASTREVIEWED.ID);

            CivicUser userLastReviewed = lastReviewed.user();
            int idLastReviewedUser = context.insertInto(CIVICLASTREVIEWEDUSER,
                    CIVICLASTREVIEWEDUSER.USERNAME,
                    CIVICLASTREVIEWEDUSER.NAME,
                    CIVICLASTREVIEWEDUSER.DISPLAYNAME,
                    CIVICLASTREVIEWEDUSER.ROLE,
                    CIVICLASTREVIEWEDUSER.AFFILIATION,
                    CIVICLASTREVIEWEDUSER.FEATUREDEXPERT,
                    CIVICLASTREVIEWEDUSER.AREAOFEXPERTISE,
                    CIVICLASTREVIEWEDUSER.BIO,
                    CIVICLASTREVIEWEDUSER.URL,
                    CIVICLASTREVIEWEDUSER.CREATEDAT,
                    CIVICLASTREVIEWEDUSER.LASTSEENAT,
                    CIVICLASTREVIEWEDUSER.AVATARURL,
                    CIVICLASTREVIEWEDUSER.TWITTERHANDLE,
                    CIVICLASTREVIEWEDUSER.FACEBOOKPROFILE,
                    CIVICLASTREVIEWEDUSER.LINKEDINPROFILE,
                    CIVICLASTREVIEWEDUSER.ORCID,
                    CIVICLASTREVIEWEDUSER.SIGNUPCOMPLETE,
                    CIVICLASTREVIEWEDUSER.ACCEPTEDLICENSE,
                    CIVICLASTREVIEWEDUSER.IDUSER,
                    CIVICLASTREVIEWEDUSER.CIVICLASTREVIEWEDID)
                    .values(userLastReviewed.username(),
                            userLastReviewed.name(),
                            userLastReviewed.displayName(),
                            userLastReviewed.role(),
                            userLastReviewed.affiliation(),
                            userLastReviewed.featuredExpert(),
                            userLastReviewed.areaOfExpertise(),
                            userLastReviewed.bio(),
                            userLastReviewed.url(),
                            userLastReviewed.createdAt(),
                            userLastReviewed.lastSeenAt(),
                            userLastReviewed.avatarUrl(),
                            userLastReviewed.twitterHandle(),
                            userLastReviewed.facebookProfile(),
                            userLastReviewed.linkedinProfile(),
                            userLastReviewed.orcid(),
                            userLastReviewed.signupComplete(),
                            userLastReviewed.acceptedLicense(),
                            userLastReviewed.id(),
                            idLastReviewed)
                    .returning(CIVICLASTREVIEWEDUSER.ID)
                    .fetchOne()
                    .getValue(CIVICLASTREVIEWEDUSER.ID);

            context.insertInto(CIVICLASTREVIEWEDAVATARS,
                    CIVICLASTREVIEWEDAVATARS.X14,
                    CIVICLASTREVIEWEDAVATARS.X32,
                    CIVICLASTREVIEWEDAVATARS.X64,
                    CIVICLASTREVIEWEDAVATARS.X128,
                    CIVICLASTREVIEWEDAVATARS.CIVICLASTREVIEWEDID)
                    .values(userLastReviewed.avatars().x14(),
                            userLastReviewed.avatars().x32(),
                            userLastReviewed.avatars().x64(),
                            userLastReviewed.avatars().x128(),
                            idLastReviewedUser)
                    .execute();

            int idLastReviewedOrganization = context.insertInto(CIVICLASTREVIEWEDORGANIZATION,
                    CIVICLASTREVIEWEDORGANIZATION.NAME,
                    CIVICLASTREVIEWEDORGANIZATION.URL,
                    CIVICLASTREVIEWEDORGANIZATION.IDORGANIZATION,
                    CIVICLASTREVIEWEDORGANIZATION.DESCRIPTION,
                    CIVICLASTREVIEWEDORGANIZATION.CIVICLASTREVIEWEDUSERID)
                    .values(userLastReviewed.organization().name(),
                            userLastReviewed.organization().url(),
                            userLastReviewed.organization().id(),
                            userLastReviewed.organization().description(),
                            idLastReviewedUser)
                    .returning(CIVICLASTREVIEWEDORGANIZATION.ID)
                    .fetchOne()
                    .getValue(CIVICLASTREVIEWEDORGANIZATION.ID);

            CivicProfileImage userLastReviewedProfileImage = userLastReviewed.organization().profileImage();
            if (userLastReviewedProfileImage != null) {
                context.insertInto(CIVICLASTREVIEWEDPROFILEIMAGE,
                        CIVICLASTREVIEWEDPROFILEIMAGE.X14,
                        CIVICLASTREVIEWEDPROFILEIMAGE.X32,
                        CIVICLASTREVIEWEDPROFILEIMAGE.X64,
                        CIVICLASTREVIEWEDPROFILEIMAGE.X128,
                        CIVICLASTREVIEWEDPROFILEIMAGE.X256,
                        CIVICLASTREVIEWEDPROFILEIMAGE.CIVICLASTREVIEWEDORGANIZATIONID)
                        .values(userLastReviewedProfileImage.x14(),
                                userLastReviewedProfileImage.x32(),
                                userLastReviewedProfileImage.x64(),
                                userLastReviewedProfileImage.x128(),
                                userLastReviewedProfileImage.x256(),
                                idLastReviewedOrganization)
                        .execute();
            }
        }
    }

    static void deleteAll(@NotNull DSLContext context) {
        // Start with deleting the life cycle action tables
        context.deleteFrom(CIVICLASTCOMMENTEDONPROFILEIMAGE).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONORGANIZATION).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONAVATARS).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDONUSER).execute();
        context.deleteFrom(CIVICLASTCOMMENTEDON).execute();
        context.deleteFrom(CIVICLASTMODIFIEDPROFILEIMAGE).execute();
        context.deleteFrom(CIVICLASTMODIFIEDORGANIZATION).execute();
        context.deleteFrom(CIVICLASTMODIFIEDAVATARS).execute();
        context.deleteFrom(CIVICLASTMODIFIEDUSER).execute();
        context.deleteFrom(CIVICLASTMODIFIED).execute();
        context.deleteFrom(CIVICLASTREVIEWEDPROFILEIMAGE).execute();
        context.deleteFrom(CIVICLASTREVIEWEDORGANIZATION).execute();
        context.deleteFrom(CIVICLASTREVIEWEDAVATARS).execute();
        context.deleteFrom(CIVICLASTREVIEWEDUSER).execute();
        context.deleteFrom(CIVICLASTREVIEWED).execute();
        context.deleteFrom(CIVICLIFECYCLEACTIONS).execute();

        // Then delete the source information
        context.deleteFrom(CIVICPUBLICATION).execute();
        context.deleteFrom(CIVICCLINICALTRIAL).execute();
        context.deleteFrom(CIVICSOURCE).execute();

        // Then delete the evidence items
        context.deleteFrom(CIVICEVIDENCEITEMCLINICALTRIAL).execute();
        context.deleteFrom(CIVICEVIDENCEITEMPUBLICATION).execute();
        context.deleteFrom(CIVICEVIDENCEITEMSOURCE).execute();
        context.deleteFrom(CIVICDISEASE).execute();
        context.deleteFrom(CIVICDRUG).execute();
        context.deleteFrom(CIVICEVIDENCEITEM).execute();

        // Then delete the variant group tables.
        context.deleteFrom(CIVICVARIANTGROUPCOORDINATES).execute();
        context.deleteFrom(CIVICVARIANTGROUPTYPE).execute();
        context.deleteFrom(CIVICVARIANTGROUPVARIANT).execute();
        context.deleteFrom(CIVICVARIANTGROUP).execute();

        // Then delete the tables directly under civic
        context.deleteFrom(CIVICVARIANTALIAS).execute();
        context.deleteFrom(CIVICVARIANTTYPE).execute();
        context.deleteFrom(CIVICPROVISIONALVALUE).execute();
        context.deleteFrom(CIVICCOORDINATES).execute();
        context.deleteFrom(CIVICASSERTION).execute();
        context.deleteFrom(CIVICHGVSEXPRESSION).execute();
        context.deleteFrom(CIVICCLINVARENTRY).execute();

        // Finally, delete the civic table
        context.deleteFrom(CIVIC).execute();
    }
}

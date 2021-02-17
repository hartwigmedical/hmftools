package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.tables.Categoryvariant.CATEGORYVARIANT;
import static com.hartwig.hmftools.ckb.database.tables.Categoryvariantpath.CATEGORYVARIANTPATH;
import static com.hartwig.hmftools.ckb.database.tables.Gene.GENE;
import static com.hartwig.hmftools.ckb.database.tables.Genedescription.GENEDESCRIPTION;
import static com.hartwig.hmftools.ckb.database.tables.Genedescriptionreference.GENEDESCRIPTIONREFERENCE;
import static com.hartwig.hmftools.ckb.database.tables.Genesynonym.GENESYNONYM;
import static com.hartwig.hmftools.ckb.database.tables.Geneterm.GENETERM;
import static com.hartwig.hmftools.ckb.database.tables.Membervariant.MEMBERVARIANT;
import static com.hartwig.hmftools.ckb.database.tables.Membervariantdescription.MEMBERVARIANTDESCRIPTION;
import static com.hartwig.hmftools.ckb.database.tables.Membervariantdescriptionreference.MEMBERVARIANTDESCRIPTIONREFERENCE;
import static com.hartwig.hmftools.ckb.database.tables.Transcriptcoordinate.TRANSCRIPTCOORDINATE;
import static com.hartwig.hmftools.ckb.database.tables.Variant.VARIANT;
import static com.hartwig.hmftools.ckb.database.tables.Variantdescription.VARIANTDESCRIPTION;
import static com.hartwig.hmftools.ckb.database.tables.Variantdescriptionreference.VARIANTDESCRIPTIONREFERENCE;

import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.variant.CategoryVariant;
import com.hartwig.hmftools.ckb.datamodel.variant.CategoryVariantPath;
import com.hartwig.hmftools.ckb.datamodel.variant.Gene;
import com.hartwig.hmftools.ckb.datamodel.variant.GeneDescription;
import com.hartwig.hmftools.ckb.datamodel.variant.MemberVariant;
import com.hartwig.hmftools.ckb.datamodel.variant.TranscriptCoordinate;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;
import com.hartwig.hmftools.ckb.datamodel.variant.VariantDescription;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;

class VariantDAO {

    @NotNull
    private final DSLContext context;

    public VariantDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    public void deleteAll() {
        // Note that deletions should go from branch to root
        context.deleteFrom(MEMBERVARIANTDESCRIPTIONREFERENCE).execute();
        context.deleteFrom(MEMBERVARIANTDESCRIPTION).execute();
        context.deleteFrom(MEMBERVARIANT).execute();

        context.deleteFrom(CATEGORYVARIANT).execute();
        context.deleteFrom(CATEGORYVARIANTPATH).execute();

        context.deleteFrom(TRANSCRIPTCOORDINATE).execute();

        context.deleteFrom(VARIANTDESCRIPTIONREFERENCE).execute();
        context.deleteFrom(VARIANTDESCRIPTION).execute();

        context.deleteFrom(GENEDESCRIPTIONREFERENCE).execute();
        context.deleteFrom(GENEDESCRIPTION).execute();
        context.deleteFrom(GENESYNONYM).execute();
        context.deleteFrom(GENETERM).execute();
        context.deleteFrom(GENE).execute();

        context.deleteFrom(VARIANT).execute();
    }

    public void write(@NotNull Variant variant, int ckbEntryId) {
        int id = context.insertInto(VARIANT,
                VARIANT.CKBENTRYID,
                VARIANT.CKBVARIANTID,
                VARIANT.FULLNAME,
                VARIANT.IMPACT,
                VARIANT.PROTEINEFFECT,
                VARIANT.TYPE,
                VARIANT.VARIANT_,
                VARIANT.CREATEDATE,
                VARIANT.UPDATEDATE)
                .values(ckbEntryId,
                        variant.id(),
                        variant.fullName(),
                        variant.impact(),
                        variant.proteinEffect(),
                        variant.type(),
                        variant.variant(),
                        Util.sqlDate(variant.createDate()),
                        Util.sqlDate(variant.updateDate()))
                .returning(VARIANT.ID)
                .fetchOne()
                .getValue(VARIANT.ID);

        writeGene(variant.gene(), id);

        for (VariantDescription variantDescription : variant.descriptions()) {
            writeVariantDescription(variantDescription, id);
        }

        writeTranscriptCoordinate(variant.referenceTranscriptCoordinate(), id, true);
        for (TranscriptCoordinate transcriptCoordinate : variant.allTranscriptCoordinates()) {
            writeTranscriptCoordinate(transcriptCoordinate, id, false);
        }

        for (CategoryVariantPath categoryVariantPath : variant.categoryVariantPaths()) {
            writeCategoryVariantPath(categoryVariantPath, id);
        }

        for (MemberVariant memberVariant : variant.memberVariants()) {
            writeMemberVariant(memberVariant, id);
        }
    }

    private void writeGene(@NotNull Gene gene, int variantId) {
        int id = context.insertInto(GENE,
                GENE.VARIANTID,
                GENE.CKBGENEID,
                GENE.GENESYMBOL,
                GENE.ENTREZID,
                GENE.CHROMOSOME,
                GENE.MAPLOCATION,
                GENE.CANONICALTRANSCRIPT,
                GENE.GENEROLE,
                GENE.CREATEDATE,
                GENE.UPDATEDATE)
                .values(variantId,
                        gene.id(),
                        gene.geneSymbol(),
                        gene.entrezId(),
                        gene.chromosome(),
                        gene.mapLocation(),
                        gene.canonicalTranscript(),
                        gene.geneRole(),
                        Util.sqlDate(gene.createDate()),
                        Util.sqlDate(gene.updateDate()))
                .returning(GENE.ID)
                .fetchOne()
                .getValue(GENE.ID);

        for (String term : gene.terms()) {
            context.insertInto(GENETERM, GENETERM.GENEID, GENETERM.TERM).values(id, term).execute();
        }

        for (String synonym : gene.synonyms()) {
            context.insertInto(GENESYNONYM, GENESYNONYM.GENEID, GENESYNONYM.SYNONYM).values(id, synonym).execute();
        }

        for (GeneDescription description : gene.descriptions()) {
            writeGeneDescription(description, id);
        }
    }

    private void writeGeneDescription(@NotNull GeneDescription geneDescription, int geneId) {
        int id = context.insertInto(GENEDESCRIPTION, GENEDESCRIPTION.GENEID, GENEDESCRIPTION.DESCRIPTION)
                .values(geneId, geneDescription.description())
                .returning(GENEDESCRIPTION.ID)
                .fetchOne()
                .getValue(GENEDESCRIPTION.ID);

        for (Reference reference : geneDescription.references()) {
            context.insertInto(GENEDESCRIPTIONREFERENCE,
                    GENEDESCRIPTIONREFERENCE.GENEDESCRIPTIONID,
                    GENEDESCRIPTIONREFERENCE.CKBREFERENCEID,
                    GENEDESCRIPTIONREFERENCE.PUBMEDID,
                    GENEDESCRIPTIONREFERENCE.TITLE,
                    GENEDESCRIPTIONREFERENCE.ABSTRACTTEXT,
                    GENEDESCRIPTIONREFERENCE.URL,
                    GENEDESCRIPTIONREFERENCE.AUTHORS,
                    GENEDESCRIPTIONREFERENCE.JOURNAL,
                    GENEDESCRIPTIONREFERENCE.VOLUME,
                    GENEDESCRIPTIONREFERENCE.ISSUE,
                    GENEDESCRIPTIONREFERENCE.DATE,
                    GENEDESCRIPTIONREFERENCE.YEAR)
                    .values(id,
                            reference.id(),
                            reference.pubMedId(),
                            reference.title(),
                            reference.abstractText(),
                            reference.url(),
                            reference.authors(),
                            reference.journal(),
                            reference.volume(),
                            reference.issue(),
                            reference.date(),
                            reference.year())
                    .execute();
        }
    }

    private void writeVariantDescription(@NotNull VariantDescription variantDescription, int variantId) {
        int id = context.insertInto(VARIANTDESCRIPTION, VARIANTDESCRIPTION.VARIANTID, VARIANTDESCRIPTION.DESCRIPTION)
                .values(variantId, variantDescription.description())
                .returning(VARIANTDESCRIPTION.ID)
                .fetchOne()
                .getValue(VARIANTDESCRIPTION.ID);

        for (Reference reference : variantDescription.references()) {
            context.insertInto(VARIANTDESCRIPTIONREFERENCE,
                    VARIANTDESCRIPTIONREFERENCE.VARIANTDESCRIPTIONID,
                    VARIANTDESCRIPTIONREFERENCE.CKBREFERENCEID,
                    VARIANTDESCRIPTIONREFERENCE.PUBMEDID,
                    VARIANTDESCRIPTIONREFERENCE.TITLE,
                    VARIANTDESCRIPTIONREFERENCE.ABSTRACTTEXT,
                    VARIANTDESCRIPTIONREFERENCE.URL,
                    VARIANTDESCRIPTIONREFERENCE.AUTHORS,
                    VARIANTDESCRIPTIONREFERENCE.JOURNAL,
                    VARIANTDESCRIPTIONREFERENCE.VOLUME,
                    VARIANTDESCRIPTIONREFERENCE.ISSUE,
                    VARIANTDESCRIPTIONREFERENCE.DATE,
                    VARIANTDESCRIPTIONREFERENCE.YEAR)
                    .values(id,
                            reference.id(),
                            reference.pubMedId(),
                            reference.title(),
                            reference.abstractText(),
                            reference.url(),
                            reference.authors(),
                            reference.journal(),
                            reference.volume(),
                            reference.issue(),
                            reference.date(),
                            reference.year())
                    .execute();
        }
    }

    private void writeTranscriptCoordinate(@Nullable TranscriptCoordinate transcriptCoordinate, int variantId,
            boolean isReferenceTranscriptCoordinate) {
        if (transcriptCoordinate != null) {
            context.insertInto(TRANSCRIPTCOORDINATE,
                    TRANSCRIPTCOORDINATE.VARIANTID,
                    TRANSCRIPTCOORDINATE.ISREFERENCETRANSCRIPTCOORDINATE,
                    TRANSCRIPTCOORDINATE.TRANSCRIPT,
                    TRANSCRIPTCOORDINATE.GDNA,
                    TRANSCRIPTCOORDINATE.CDNA,
                    TRANSCRIPTCOORDINATE.PROTEIN,
                    TRANSCRIPTCOORDINATE.SOURCEDB,
                    TRANSCRIPTCOORDINATE.REFGENOMEBUILD)
                    .values(variantId,
                            Util.toByte(isReferenceTranscriptCoordinate),
                            transcriptCoordinate.transcript(),
                            transcriptCoordinate.gDna(),
                            transcriptCoordinate.cDna(),
                            transcriptCoordinate.protein(),
                            transcriptCoordinate.sourceDb(),
                            transcriptCoordinate.refGenomeBuild())
                    .execute();
        }
    }

    private void writeCategoryVariantPath(@NotNull CategoryVariantPath categoryVariantPath, int variantId) {
        int id = context.insertInto(CATEGORYVARIANTPATH, CATEGORYVARIANTPATH.VARIANTID, CATEGORYVARIANTPATH.VARIANTPATH)
                .values(variantId, categoryVariantPath.variantPath())
                .returning(CATEGORYVARIANTPATH.ID)
                .fetchOne()
                .getValue(CATEGORYVARIANTPATH.ID);

        for (CategoryVariant categoryVariant : categoryVariantPath.variants()) {
            context.insertInto(CATEGORYVARIANT,
                    CATEGORYVARIANT.CATEGORYVARIANTPATHID,
                    CATEGORYVARIANT.CKBVARIANTID,
                    CATEGORYVARIANT.FULLNAME,
                    CATEGORYVARIANT.IMPACT,
                    CATEGORYVARIANT.PROTEINEFFECT)
                    .values(id, categoryVariant.id(), categoryVariant.fullName(), categoryVariant.impact(), categoryVariant.proteinEffect())
                    .execute();
        }
    }

    private void writeMemberVariant(@NotNull MemberVariant memberVariant, int variantId) {
        int id = context.insertInto(MEMBERVARIANT,
                MEMBERVARIANT.VARIANTID,
                MEMBERVARIANT.CKBVARIANTID,
                MEMBERVARIANT.FULLNAME,
                MEMBERVARIANT.IMPACT,
                MEMBERVARIANT.PROTEINEFFECT)
                .values(variantId, memberVariant.id(), memberVariant.fullName(), memberVariant.impact(), memberVariant.proteinEffect())
                .returning(MEMBERVARIANT.ID)
                .fetchOne()
                .getValue(MEMBERVARIANT.ID);

        for (VariantDescription memberVariantDescription : memberVariant.descriptions()) {
            writeMemberVariantDescription(memberVariantDescription, id);
        }

    }

    private void writeMemberVariantDescription(@NotNull VariantDescription memberVariantDescription, int memberVariantId) {
        int id =
                context.insertInto(MEMBERVARIANTDESCRIPTION, MEMBERVARIANTDESCRIPTION.MEMBERVARIANTID, MEMBERVARIANTDESCRIPTION.DESCRIPTION)
                        .values(memberVariantId, memberVariantDescription.description())
                        .returning(MEMBERVARIANTDESCRIPTION.ID)
                        .fetchOne()
                        .getValue(MEMBERVARIANTDESCRIPTION.ID);

        for (Reference reference : memberVariantDescription.references()) {
            context.insertInto(MEMBERVARIANTDESCRIPTIONREFERENCE,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.MEMBERVARIANTDESCRIPTIONID,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.CKBREFERENCEID,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.PUBMEDID,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.TITLE,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.ABSTRACTTEXT,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.URL,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.AUTHORS,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.JOURNAL,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.VOLUME,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.ISSUE,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.DATE,
                    MEMBERVARIANTDESCRIPTIONREFERENCE.YEAR)
                    .values(id,
                            reference.id(),
                            reference.pubMedId(),
                            reference.title(),
                            reference.abstractText(),
                            reference.url(),
                            reference.authors(),
                            reference.journal(),
                            reference.volume(),
                            reference.issue(),
                            reference.date(),
                            reference.year())
                    .execute();
        }
    }
}

package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.tables.Categoryvariantpath.CATEGORYVARIANTPATH;
import static com.hartwig.hmftools.ckb.database.tables.Gene.GENE;
import static com.hartwig.hmftools.ckb.database.tables.Genedescription.GENEDESCRIPTION;
import static com.hartwig.hmftools.ckb.database.tables.Genedescriptionreference.GENEDESCRIPTIONREFERENCE;
import static com.hartwig.hmftools.ckb.database.tables.Genesynonym.GENESYNONYM;
import static com.hartwig.hmftools.ckb.database.tables.Geneterm.GENETERM;
import static com.hartwig.hmftools.ckb.database.tables.Membervariant.MEMBERVARIANT;
import static com.hartwig.hmftools.ckb.database.tables.Transcriptcoordinate.TRANSCRIPTCOORDINATE;
import static com.hartwig.hmftools.ckb.database.tables.Variant.VARIANT;
import static com.hartwig.hmftools.ckb.database.tables.Variantdescription.VARIANTDESCRIPTION;
import static com.hartwig.hmftools.ckb.database.tables.Variantdescriptionreference.VARIANTDESCRIPTIONREFERENCE;

import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
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
        context.deleteFrom(MEMBERVARIANT).execute();
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
                VARIANT.CREATEDATE,
                VARIANT.UPDATEDATE,
                VARIANT.FULLNAME,
                VARIANT.VARIANT_,
                VARIANT.IMPACT,
                VARIANT.PROTEINEFFECT,
                VARIANT.TYPE)
                .values(ckbEntryId,
                        variant.id(),
                        Util.sqlDate(variant.createDate()),
                        Util.sqlDate(variant.updateDate()),
                        variant.fullName(),
                        variant.variant(),
                        variant.impact(),
                        variant.proteinEffect(),
                        variant.type())
                .returning(VARIANT.ID)
                .fetchOne()
                .getValue(VARIANT.ID);

        writeGene(variant.gene(), id);

        writeTranscriptCoordinate(variant.referenceTranscriptCoordinate(), id, true);
        for (TranscriptCoordinate transcriptCoordinate : variant.allTranscriptCoordinates()) {
            writeTranscriptCoordinate(transcriptCoordinate, id, false);
        }

        for (String categoryVariantPath : variant.categoryVariantPaths()) {
            context.insertInto(CATEGORYVARIANTPATH, CATEGORYVARIANTPATH.VARIANTID, CATEGORYVARIANTPATH.VARIANTPATH)
                    .values(id, categoryVariantPath)
                    .execute();
        }

        for (MemberVariant memberVariant : variant.memberVariants()) {
            writeMemberVariant(memberVariant, id);
        }

        for (VariantDescription variantDescription : variant.descriptions()) {
            writeVariantDescription(variantDescription, id);
        }
    }

    private void writeGene(@NotNull Gene gene, int variantId) {
        int id = context.insertInto(GENE,
                GENE.VARIANTID,
                GENE.CKBGENEID,
                GENE.CREATEDATE,
                GENE.UPDATEDATE,
                GENE.GENESYMBOL,
                GENE.GENEROLE,
                GENE.ENTREZID,
                GENE.CHROMOSOME,
                GENE.MAPLOCATION,
                GENE.CANONICALTRANSCRIPT)
                .values(variantId,
                        gene.id(),
                        Util.sqlDate(gene.createDate()),
                        Util.sqlDate(gene.updateDate()),
                        gene.geneSymbol(),
                        gene.geneRole(),
                        gene.entrezId(),
                        gene.chromosome(),
                        gene.mapLocation(),
                        gene.canonicalTranscript())
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
                    GENEDESCRIPTIONREFERENCE.JOURNAL,
                    GENEDESCRIPTIONREFERENCE.AUTHORS,
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
                            reference.journal(),
                            reference.authors(),
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
                    VARIANTDESCRIPTIONREFERENCE.JOURNAL,
                    VARIANTDESCRIPTIONREFERENCE.AUTHORS,
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
                            reference.journal(),
                            reference.authors(),
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

    private void writeMemberVariant(@NotNull MemberVariant memberVariant, int variantId) {
        context.insertInto(MEMBERVARIANT,
                MEMBERVARIANT.VARIANTID,
                MEMBERVARIANT.CKBVARIANTID,
                MEMBERVARIANT.FULLNAME,
                MEMBERVARIANT.IMPACT,
                MEMBERVARIANT.PROTEINEFFECT)
                .values(variantId, memberVariant.id(), memberVariant.fullName(), memberVariant.impact(), memberVariant.proteinEffect())
                .execute();
    }
}

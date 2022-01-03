package com.hartwig.hmftools.ckb.dao;

import static com.hartwig.hmftools.ckb.database.Tables.GENEREFERENCE;
import static com.hartwig.hmftools.ckb.database.Tables.VARIANTREFERENCE;
import static com.hartwig.hmftools.ckb.database.tables.Categoryvariantpath.CATEGORYVARIANTPATH;
import static com.hartwig.hmftools.ckb.database.tables.Gene.GENE;
import static com.hartwig.hmftools.ckb.database.tables.Genesynonym.GENESYNONYM;
import static com.hartwig.hmftools.ckb.database.tables.Geneterm.GENETERM;
import static com.hartwig.hmftools.ckb.database.tables.Membervariant.MEMBERVARIANT;
import static com.hartwig.hmftools.ckb.database.tables.Transcriptcoordinate.TRANSCRIPTCOORDINATE;
import static com.hartwig.hmftools.ckb.database.tables.Variant.VARIANT;

import com.hartwig.hmftools.ckb.datamodel.reference.Reference;
import com.hartwig.hmftools.ckb.datamodel.variant.Gene;
import com.hartwig.hmftools.ckb.datamodel.variant.MemberVariant;
import com.hartwig.hmftools.ckb.datamodel.variant.TranscriptCoordinate;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

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
        context.deleteFrom(VARIANTREFERENCE).execute();

        context.deleteFrom(GENEREFERENCE).execute();
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
                VARIANT.TYPE,
                VARIANT.DESCRIPTION)
                .values(ckbEntryId,
                        variant.id(),
                        variant.createDate(),
                        variant.updateDate(),
                        variant.fullName(),
                        variant.variant(),
                        variant.impact(),
                        variant.proteinEffect(),
                        variant.type(),
                        variant.description())
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

        for (Reference variantReference : variant.references()) {
            writeVariantReference(variantReference, id);
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
                GENE.CANONICALTRANSCRIPT,
                GENE.DESCRIPTION)
                .values(variantId,
                        gene.id(),
                        gene.createDate(),
                        gene.updateDate(),
                        gene.geneSymbol(),
                        gene.geneRole(),
                        gene.entrezId(),
                        gene.chromosome(),
                        gene.mapLocation(),
                        gene.canonicalTranscript(),
                        gene.description())
                .returning(GENE.ID)
                .fetchOne()
                .getValue(GENE.ID);

        for (String term : gene.terms()) {
            context.insertInto(GENETERM, GENETERM.GENEID, GENETERM.TERM).values(id, term).execute();
        }

        for (String synonym : gene.synonyms()) {
            context.insertInto(GENESYNONYM, GENESYNONYM.GENEID, GENESYNONYM.SYNONYM).values(id, synonym).execute();
        }

        for (Reference geneReference : gene.references()) {
            writeGeneReference(geneReference, id);
        }
    }

    private void writeGeneReference(@NotNull Reference reference, int geneId) {
        context.insertInto(GENEREFERENCE,
                GENEREFERENCE.GENEID,
                GENEREFERENCE.CKBREFERENCEID,
                GENEREFERENCE.PUBMEDID,
                GENEREFERENCE.TITLE,
                GENEREFERENCE.ABSTRACTTEXT,
                GENEREFERENCE.URL,
                GENEREFERENCE.JOURNAL,
                GENEREFERENCE.AUTHORS,
                GENEREFERENCE.VOLUME,
                GENEREFERENCE.ISSUE,
                GENEREFERENCE.DATE,
                GENEREFERENCE.YEAR)
                .values(geneId,
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

    private void writeVariantReference(@NotNull Reference reference, int variantId) {
        context.insertInto(VARIANTREFERENCE,
                VARIANTREFERENCE.VARIANTID,
                VARIANTREFERENCE.CKBREFERENCEID,
                VARIANTREFERENCE.PUBMEDID,
                VARIANTREFERENCE.TITLE,
                VARIANTREFERENCE.ABSTRACTTEXT,
                VARIANTREFERENCE.URL,
                VARIANTREFERENCE.JOURNAL,
                VARIANTREFERENCE.AUTHORS,
                VARIANTREFERENCE.VOLUME,
                VARIANTREFERENCE.ISSUE,
                VARIANTREFERENCE.DATE,
                VARIANTREFERENCE.YEAR)
                .values(variantId,
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

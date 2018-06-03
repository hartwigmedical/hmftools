package com.hartwig.hmftools.svannotation;

import static org.ensembl.database.homo_sapiens_core.Tables.COORD_SYSTEM;
import static org.ensembl.database.homo_sapiens_core.Tables.EXON;
import static org.ensembl.database.homo_sapiens_core.Tables.EXON_TRANSCRIPT;
import static org.ensembl.database.homo_sapiens_core.Tables.GENE;
import static org.ensembl.database.homo_sapiens_core.Tables.KARYOTYPE;
import static org.ensembl.database.homo_sapiens_core.Tables.OBJECT_XREF;
import static org.ensembl.database.homo_sapiens_core.Tables.SEQ_REGION;
import static org.ensembl.database.homo_sapiens_core.Tables.TRANSCRIPT;
import static org.ensembl.database.homo_sapiens_core.Tables.XREF;
import static org.jooq.impl.DSL.decode;
import static org.jooq.impl.DSL.groupConcatDistinct;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.svannotation.annotations.GeneAnnotation;
import com.hartwig.hmftools.svannotation.annotations.StructuralVariantAnnotation;
import com.hartwig.hmftools.svannotation.annotations.Transcript;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ensembl.database.homo_sapiens_core.enums.GeneStatus;
import org.ensembl.database.homo_sapiens_core.enums.ObjectXrefEnsemblObjectType;
import org.ensembl.database.homo_sapiens_core.tables.Xref;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.Condition;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;
import org.jooq.types.UInteger;

public class MySQLAnnotator implements VariantAnnotator {

    private static final Logger LOGGER = LogManager.getLogger(MySQLAnnotator.class);

    private static final String ENTREZ_IDS = "ENTREZ_IDS";
    private static final String KARYOTYPE_BAND = "KARYOTYPE_BAND";

    @NotNull
    private final DSLContext context;
    @NotNull
    private final UInteger coordSystemId;

    @NotNull
    public static VariantAnnotator make(@NotNull String url) throws SQLException {
        return new MySQLAnnotator(url);
    }

    private MySQLAnnotator(@NotNull String url) throws SQLException {
        System.setProperty("org.jooq.no-logo", "true");
        final Connection conn = DriverManager.getConnection(url);
        context = DSL.using(conn, SQLDialect.MYSQL);
        coordSystemId = findCoordSystemId();
    }

    public MySQLAnnotator(@NotNull final DSLContext dbConnection)
    {
        context = dbConnection;
        coordSystemId = findCoordSystemId();
    }

    @NotNull
    private UInteger findCoordSystemId() {
        return context.select(COORD_SYSTEM.COORD_SYSTEM_ID)
                .from(COORD_SYSTEM)
                .where(COORD_SYSTEM.VERSION.eq("GRCh37"))
                .orderBy(COORD_SYSTEM.RANK)
                .limit(1)
                .fetchOne()
                .value1();
    }

    @Override
    @NotNull
    public List<StructuralVariantAnnotation> annotateVariants(@NotNull List<EnrichedStructuralVariant> variants) {

        List<StructuralVariantAnnotation> annotatedVars = Lists.newArrayList();

//        for(final EnrichedStructuralVariant variant : variants)
//        {
//             annotatedVars.add(annotateVariant(variant));
//        }
//
//        return annotatedVars;
        return variants.stream().map(this::annotateVariant).collect(Collectors.toList());
    }

    @NotNull
    private StructuralVariantAnnotation annotateVariant(@NotNull EnrichedStructuralVariant variant) {
        final StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(variant);
        annotation.annotations().addAll(annotateBreakend(variant, true, variant.start().chromosome(), variant.start().position()));
        annotation.annotations().addAll(annotateBreakend(variant, false, variant.end().chromosome(), variant.end().position()));
        return annotation;
    }

    @NotNull
    private List<GeneAnnotation> annotateBreakend(@NotNull EnrichedStructuralVariant variant, final boolean isStart,
            @NotNull String chromosome, final long position) {
        final List<GeneAnnotation> result = Lists.newArrayList();

        final Result<?> genes = queryGenesOnChromosomeAndPosition(chromosome, position);

        for (final Record gene : genes) {
            final UInteger geneId = gene.get(GENE.GENE_ID);
            final String geneName = gene.get(XREF.DISPLAY_LABEL);
            final String geneStableId = gene.get(GENE.STABLE_ID);
            final UInteger canonicalTranscriptId = gene.get(GENE.CANONICAL_TRANSCRIPT_ID);
            final int geneStrand = gene.get(GENE.SEQ_REGION_STRAND);

            final String entrezIdsStr = gene.get(ENTREZ_IDS, String.class);

//            if(entrezIdsStr == null) {
//                LOGGER.debug("var({}) missing an entrezId", variant.id());
//            }

            final List<Integer> entrezIds = (entrezIdsStr == null || entrezIdsStr.isEmpty()) ? Lists.newArrayList() :
                Arrays.stream(entrezIdsStr.split(",")).map(Integer::parseInt).collect(Collectors.toList());

            final String karyotypeBand = gene.get(KARYOTYPE_BAND, String.class);

            final List<String> synonyms = context.select(XREF.DBPRIMARY_ACC)
                    .from(XREF)
                    .innerJoin(OBJECT_XREF)
                    .on(OBJECT_XREF.XREF_ID.eq(XREF.XREF_ID))
                    .and(OBJECT_XREF.ENSEMBL_ID.eq(geneId))
                    .and(OBJECT_XREF.ENSEMBL_OBJECT_TYPE.eq(ObjectXrefEnsemblObjectType.Gene))
                    .fetch()
                    .stream()
                    .map(r -> r.get(XREF.DBPRIMARY_ACC))
                    .collect(Collectors.toList());

            final GeneAnnotation geneAnnotation =
                    new GeneAnnotation(variant, isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand);

            final Result<?> transcripts = context.select(TRANSCRIPT.TRANSCRIPT_ID, TRANSCRIPT.STABLE_ID)
                    .from(TRANSCRIPT)
                    .where(TRANSCRIPT.GENE_ID.eq(geneId))
                    .fetch();

            for (final Record transcriptRecord : transcripts) {
                Transcript transcript = buildTranscript(geneAnnotation, transcriptRecord, position, canonicalTranscriptId, geneStrand > 0);

                if (transcript != null) {
                    geneAnnotation.addTranscript(transcript);
                }
            }

            if (!geneAnnotation.transcripts().isEmpty()) {
                result.add(geneAnnotation);
            }
        }

        return result;
    }

    @NotNull
    private Result<?> queryGenesOnChromosomeAndPosition(@NotNull String chromosome, long position) {
        final int promoterDistance = 10000;
        final byte zero = 0;
        final Xref ENTREZ_XREF = XREF.as("entrez_xref");
        return context.select(GENE.GENE_ID,
                XREF.DISPLAY_LABEL,
                GENE.STABLE_ID,
                GENE.CANONICAL_TRANSCRIPT_ID,
                GENE.SEQ_REGION_STRAND,
                groupConcatDistinct(ENTREZ_XREF.DBPRIMARY_ACC).separator(",").as(ENTREZ_IDS),
                groupConcatDistinct(KARYOTYPE.BAND).separator("-").as(KARYOTYPE_BAND))
                .from(GENE)
                .innerJoin(SEQ_REGION)
                .on(GENE.SEQ_REGION_ID.eq(SEQ_REGION.SEQ_REGION_ID))
                .and(SEQ_REGION.NAME.eq(chromosome))
                .and(SEQ_REGION.COORD_SYSTEM_ID.eq(coordSystemId))
                .innerJoin(KARYOTYPE)
                .on(GENE.SEQ_REGION_ID.eq(KARYOTYPE.SEQ_REGION_ID))
                .innerJoin(OBJECT_XREF)
                .on(GENE.GENE_ID.eq(OBJECT_XREF.ENSEMBL_ID))
                .and(OBJECT_XREF.ENSEMBL_OBJECT_TYPE.eq(ObjectXrefEnsemblObjectType.Gene))
                .leftJoin(ENTREZ_XREF) // was an inner join before
                .on(OBJECT_XREF.XREF_ID.eq(ENTREZ_XREF.XREF_ID))
                .and(ENTREZ_XREF.EXTERNAL_DB_ID.eq(UInteger.valueOf(1300)))
                .innerJoin(XREF)
                .on(XREF.XREF_ID.eq(GENE.DISPLAY_XREF_ID))
                .where((GENE.STATUS.eq(GeneStatus.KNOWN).or(GENE.STATUS.eq(GeneStatus.NOVEL))))
                .and(decode().when(GENE.SEQ_REGION_STRAND.gt(zero),
                        decode().when(GENE.SEQ_REGION_START.ge(UInteger.valueOf(promoterDistance)),
                                GENE.SEQ_REGION_START.sub(promoterDistance)).otherwise(GENE.SEQ_REGION_START))
                        .otherwise(GENE.SEQ_REGION_START)
                        .le(UInteger.valueOf(position)))
                .and(decode().when(GENE.SEQ_REGION_STRAND.lt(zero), GENE.SEQ_REGION_END.add(promoterDistance))
                        .otherwise(GENE.SEQ_REGION_END)
                        .ge(UInteger.valueOf(position)))
                .and(geneStartInKaryotypeBand().or(geneEndInKaryotypeBand()))
                .groupBy(GENE.GENE_ID)
                .fetch();
    }

    @Nullable
    private Transcript buildTranscript(@NotNull GeneAnnotation parent, @NotNull Record transcript, long position,
            @NotNull UInteger canonicalTranscriptId, boolean isForwardStrand) {
        final UInteger transcriptId = transcript.get(TRANSCRIPT.TRANSCRIPT_ID);
        final boolean canonical = transcriptId.equals(canonicalTranscriptId);
        final String transcriptStableId = transcript.get(TRANSCRIPT.STABLE_ID);

        final Record exonLeft = context.select(EXON_TRANSCRIPT.RANK, EXON.PHASE, EXON.END_PHASE)
                .from(EXON_TRANSCRIPT)
                .innerJoin(EXON)
                .on(EXON.EXON_ID.eq(EXON_TRANSCRIPT.EXON_ID))
                .where(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(transcriptId))
                .and(EXON.SEQ_REGION_START.le(UInteger.valueOf(position)))
                .orderBy(EXON.SEQ_REGION_START.desc())
                .limit(1)
                .fetchOne();

        final Record exonRight = context.select(EXON_TRANSCRIPT.RANK, EXON.PHASE, EXON.END_PHASE)
                .from(EXON_TRANSCRIPT)
                .innerJoin(EXON)
                .on(EXON.EXON_ID.eq(EXON_TRANSCRIPT.EXON_ID))
                .where(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(transcriptId))
                .and(EXON.SEQ_REGION_END.ge(UInteger.valueOf(position)))
                .orderBy(EXON.SEQ_REGION_END.asc())
                .limit(1)
                .fetchOne();

        final int exonMax = context.select(EXON_TRANSCRIPT.RANK)
                .from(EXON_TRANSCRIPT)
                .where(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(transcriptId))
                .orderBy(EXON_TRANSCRIPT.RANK.desc())
                .limit(1)
                .fetchOne()
                .value1();

        final int exonUpstream;
        final int exonUpstreamPhase;
        final int exonDownstream;
        final int exonDownstreamPhase;

        if (isForwardStrand) {
            exonUpstream = exonLeft == null ? 0 : exonLeft.get(EXON_TRANSCRIPT.RANK);
            exonUpstreamPhase = exonLeft == null ? -1 : exonLeft.get(EXON.END_PHASE);
            exonDownstream = exonRight == null ? 0 : exonRight.get(EXON_TRANSCRIPT.RANK);
            exonDownstreamPhase = exonRight == null ? -1 : exonRight.get(EXON.PHASE);
        } else {
            exonDownstream = exonLeft == null ? 0 : exonLeft.get(EXON_TRANSCRIPT.RANK);
            exonDownstreamPhase = exonLeft == null ? -1 : exonLeft.get(EXON.PHASE);
            exonUpstream = exonRight == null ? 0 : exonRight.get(EXON_TRANSCRIPT.RANK);
            exonUpstreamPhase = exonRight == null ? -1 : exonRight.get(EXON.END_PHASE);
        }

        if (exonUpstream > 0 && exonDownstream == 0) {
            // NERA: past the last exon
            return null;
        } else {
            return new Transcript(parent,
                    transcriptStableId,
                    exonUpstream,
                    exonUpstreamPhase,
                    exonDownstream,
                    exonDownstreamPhase,
                    exonMax,
                    canonical);
        }
    }

    @NotNull
    private static Condition geneStartInKaryotypeBand() {
        return GENE.SEQ_REGION_START.ge(KARYOTYPE.SEQ_REGION_START).and(GENE.SEQ_REGION_START.le(KARYOTYPE.SEQ_REGION_END));
    }

    @NotNull
    private static Condition geneEndInKaryotypeBand() {
        return GENE.SEQ_REGION_END.ge(KARYOTYPE.SEQ_REGION_START).and(GENE.SEQ_REGION_END.le(KARYOTYPE.SEQ_REGION_END));
    }
}

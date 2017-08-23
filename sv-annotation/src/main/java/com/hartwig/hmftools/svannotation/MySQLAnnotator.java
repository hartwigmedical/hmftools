package com.hartwig.hmftools.svannotation;

import static org.jooq.impl.DSL.decode;
import static org.jooq.impl.DSL.field;
import static org.jooq.impl.DSL.name;
import static org.jooq.impl.DSL.table;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.structural.StructuralVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;

public class MySQLAnnotator implements StructuralVariantAnnotator {

    private static final Logger LOGGER = LogManager.getLogger(MySQLAnnotator.class);

    private final DSLContext context;
    private final int coord_system_id;

    public static StructuralVariantAnnotator make(final String url) throws SQLException {
        return new MySQLAnnotator(url);
    }

    private MySQLAnnotator(final String url) throws SQLException {
        System.setProperty("org.jooq.no-logo", "true");
        final Connection conn = DriverManager.getConnection(url);
        context = DSL.using(conn, SQLDialect.MYSQL);
        coord_system_id = findCoordSystemId();
    }

    private int findCoordSystemId() {
        return context.select(field(name("coord_system_id")))
                .from(table(name("coord_system")))
                .where(field(name("version")).eq("GRCh37"))
                .orderBy(field(name("rank")))
                .limit(1)
                .fetchOne(0, Integer.class);
    }

    @Override
    public List<StructuralVariantAnnotation> annotateVariants(final List<StructuralVariant> variants) {
        return variants.stream().map(this::annotateVariant).collect(Collectors.toList());
    }

    private StructuralVariantAnnotation annotateVariant(final StructuralVariant variant) {
        final StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(variant);

        annotation.setBreakendAnnotations(
                annotateBreakend(annotation, variant.startChromosome(), variant.startPosition(), variant.startOrientation(),
                        variant.startAF()),
                annotateBreakend(annotation, variant.endChromosome(), variant.endPosition(), variant.endOrientation(), variant.endAF()));

        return annotation;
    }

    private BreakendAnnotations annotateBreakend(final StructuralVariantAnnotation parent, String chromosome, final long position,
            final int orientation, final Double alleleFrequency) {
        final BreakendAnnotations annotation = new BreakendAnnotations(parent, chromosome, position, orientation, alleleFrequency);

        final int PROMOTER_DISTANCE = 10000;

        // start with the overlapping genes
        final Result<?> genes =
                context.select(field(name("gene", "gene_id")), field(name("xref", "display_label")), field(name("gene", "stable_id")),
                        field(name("gene", "canonical_transcript_id")), field(name("gene", "seq_region_strand")))
                        .from(table(name("gene")))
                        .innerJoin(table(name("seq_region")))
                        .on(field(name("gene", "seq_region_id")).eq(field(name("seq_region", "seq_region_id"))))
                        .and(field(name("seq_region", "name")).eq(chromosome))
                        .and(field(name("seq_region", "coord_system_id")).eq(coord_system_id))
                        .innerJoin(table(name("xref")))
                        .on(field(name("xref", "xref_id")).eq(field(name("gene", "display_xref_id"))))
                        .where(field(name("gene", "status")).eq("KNOWN"))
                        .and(decode().when(field(name("gene", "seq_region_strand")).ge(0),
                                field(name("gene", "seq_region_start")).sub(PROMOTER_DISTANCE))
                                .otherwise(field(name("gene", "seq_region_start")))
                                .le(position))
                        .and(decode().when(field(name("gene", "seq_region_strand")).le(0),
                                field(name("gene", "seq_region_end")).add(PROMOTER_DISTANCE))
                                .otherwise(field(name("gene", "seq_region_end")))
                                .ge(position))
                        .fetch();

        for (final Record g : genes) {
            final int gene_id = g.get(0, Integer.class);
            final String gene_name = g.get(1, String.class);
            final String gene_stable_id = g.get(2, String.class);
            final int canonical_transcript_id = g.get(3, Integer.class);
            final int gene_strand = g.get(4, Integer.class);

            final GeneAnnotation gene = new GeneAnnotation(annotation, gene_name, gene_stable_id, gene_strand);
            annotation.addGeneAnnotation(gene);

            final Result<?> transcripts = context.select(field(name("transcript_id")), field(name("stable_id")))
                    .from(table(name("transcript")))
                    .where(field(name("gene_id")).eq(gene_id))
                    .fetch();

            for (final Record t : transcripts) {
                final int transcript_id = t.get(0, Integer.class);
                final boolean canonical = transcript_id == canonical_transcript_id;
                final String transcript_stable_id = t.get(1, String.class);

                final Record exonLeft = context.select(field(name("exon_transcript", "rank")), field(name("exon", "phase")),
                        field(name("exon", "end_phase")))
                        .from(table(name("exon_transcript")))
                        .innerJoin(table(name("exon")))
                        .on(field(name("exon", "exon_id")).eq(field(name("exon_transcript", "exon_id"))))
                        .where(field(name("exon_transcript", "transcript_id")).eq(transcript_id))
                        .and(field(name("exon", "seq_region_start")).le(position))
                        .orderBy(field(name("exon", "seq_region_start")).desc())
                        .limit(1)
                        .fetchOne();

                final Record exonRight = context.select(field(name("exon_transcript", "rank")), field(name("exon", "phase")),
                        field(name("exon", "end_phase")))
                        .from(table(name("exon_transcript")))
                        .innerJoin(table(name("exon")))
                        .on(field(name("exon", "exon_id")).eq(field(name("exon_transcript", "exon_id"))))
                        .where(field(name("exon_transcript", "transcript_id")).eq(transcript_id))
                        .and(field(name("exon", "seq_region_end")).ge(position))
                        .orderBy(field(name("exon", "seq_region_end")).asc())
                        .limit(1)
                        .fetchOne();

                final int exon_max = context.select(field(name("rank")))
                        .from(table(name("exon_transcript")))
                        .where(field(name("transcript_id")).eq(transcript_id))
                        .orderBy(field(name("rank")).desc())
                        .limit(1)
                        .fetchOne(0, Integer.class);

                final int exon_upstream;
                final int exon_upstream_phase;
                final int exon_downstream;
                final int exon_downstream_phase;

                if (gene_strand > 0) {
                    // forward strand
                    exon_upstream = exonLeft == null ? 0 : exonLeft.get(0, Integer.class);
                    exon_upstream_phase = exonLeft == null ? 0 : exonLeft.get(2, Integer.class);
                    exon_downstream = exonRight == null ? 0 : exonRight.get(0, Integer.class);
                    exon_downstream_phase = exonRight == null ? 0 : exonRight.get(1, Integer.class);
                } else {
                    // reverse strand
                    exon_downstream = exonLeft == null ? 0 : exonLeft.get(0, Integer.class);
                    exon_downstream_phase = exonLeft == null ? 0 : exonLeft.get(1, Integer.class);
                    exon_upstream = exonRight == null ? 0 : exonRight.get(0, Integer.class);
                    exon_upstream_phase = exonRight == null ? 0 : exonRight.get(2, Integer.class);
                }

                final TranscriptAnnotation transcript =
                        new TranscriptAnnotation(gene, transcript_stable_id, exon_upstream, exon_upstream_phase, exon_downstream,
                                exon_downstream_phase, exon_max, canonical);
                gene.addTranscriptAnnotation(transcript);
            }
        }

        return annotation;
    }
}

package com.hartwig.hmftools.svannotation;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.variant.structural.annotation.SvGeneTranscriptCollection.PRE_GENE_PROMOTOR_DISTANCE;

import static org.ensembl.database.homo_sapiens_core.tables.CoordSystem.COORD_SYSTEM;
import static org.ensembl.database.homo_sapiens_core.tables.Exon.EXON;
import static org.ensembl.database.homo_sapiens_core.tables.ExonTranscript.EXON_TRANSCRIPT;
import static org.ensembl.database.homo_sapiens_core.tables.Gene.GENE;
import static org.ensembl.database.homo_sapiens_core.tables.Karyotype.KARYOTYPE;
import static org.ensembl.database.homo_sapiens_core.tables.ObjectXref.OBJECT_XREF;
import static org.ensembl.database.homo_sapiens_core.tables.SeqRegion.SEQ_REGION;
import static org.ensembl.database.homo_sapiens_core.tables.Transcript.TRANSCRIPT;
import static org.ensembl.database.homo_sapiens_core.tables.Xref.XREF;
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
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariantLeg;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ensembl.database.homo_sapiens_core.enums.ObjectXrefEnsemblObjectType;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;
import org.jooq.types.UInteger;

public class MySQLAnnotator implements VariantAnnotator
{
    private EnsemblDAO mEnsemblDAO;

    private static final String ENTREZ_IDS = "ENTREZ_IDS";
    private static final String KARYOTYPE_BAND = "KARYOTYPE_BAND";

    // @NotNull
    private final DSLContext mDslContext;

    private static final Logger LOGGER = LogManager.getLogger(MySQLAnnotator.class);

    @NotNull
    public static VariantAnnotator make(@NotNull String url) throws SQLException
    {
        return new MySQLAnnotator(url);
    }

    private MySQLAnnotator(@NotNull String url) throws SQLException
    {
        System.setProperty("org.jooq.no-logo", "true");
        final Connection conn = DriverManager.getConnection(url);
        final DSLContext context = DSL.using(conn, SQLDialect.MYSQL);
        mDslContext = context;
        mEnsemblDAO = new EnsemblDAO(context);
    }

    public MySQLAnnotator(@NotNull final DSLContext context)
    {
        mDslContext = context;
        mEnsemblDAO = new EnsemblDAO(context);
    }

    public final EnsemblDAO getEnsemblDAO() { return mEnsemblDAO; }

    @Override
    @NotNull
    public List<StructuralVariantAnnotation> annotateVariants(@NotNull List<EnrichedStructuralVariant> variants)
    {
        LOGGER.debug("annotating {} variants with gene and transcript info from ensembl db", variants.size());
        return variants.stream().map(this::annotateVariant).collect(Collectors.toList());
    }

    @NotNull
    private StructuralVariantAnnotation annotateVariant(@NotNull EnrichedStructuralVariant variant)
    {
        final StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(variant);

        annotation.annotations().addAll(annotateBreakend(variant, true, variant.start().chromosome(), variant.start().position()));

        EnrichedStructuralVariantLeg endLeg = variant.end();
        if (endLeg != null)
        {
            annotation.annotations().addAll(annotateBreakend(variant, false, endLeg.chromosome(), endLeg.position()));
        }

        return annotation;
    }

    @NotNull
    private List<GeneAnnotation> annotateBreakend(@NotNull EnrichedStructuralVariant variant, final boolean isStart, @NotNull String chromosome, final long position)
    {
        final List<GeneAnnotation> result = Lists.newArrayList();

        final Result<?> genes = mEnsemblDAO.queryGenesOnChromosomeAndPosition(chromosome, position, PRE_GENE_PROMOTOR_DISTANCE);

        for (final Record gene : genes)
        {
            final UInteger geneId = gene.get(GENE.GENE_ID);
            final String geneName = gene.get(XREF.DISPLAY_LABEL);
            final String geneStableId = gene.get(GENE.STABLE_ID);
            final UInteger canonicalTranscriptId = gene.get(GENE.CANONICAL_TRANSCRIPT_ID);
            final int geneStrand = gene.get(GENE.SEQ_REGION_STRAND);

            final String entrezIdsStr = gene.get(ENTREZ_IDS, String.class);

            final List<Integer> entrezIds = (entrezIdsStr == null || entrezIdsStr.isEmpty())
                    ? Lists.newArrayList()
                    : Arrays.stream(entrezIdsStr.split(",")).map(Integer::parseInt).collect(Collectors.toList());

            final String karyotypeBand = gene.get(KARYOTYPE_BAND, String.class);

            final List<String> synonyms = mDslContext.select(XREF.DBPRIMARY_ACC)
                    .from(XREF)
                    .innerJoin(OBJECT_XREF)
                    .on(OBJECT_XREF.XREF_ID.eq(XREF.XREF_ID))
                    .and(OBJECT_XREF.ENSEMBL_ID.eq(geneId))
                    .and(OBJECT_XREF.ENSEMBL_OBJECT_TYPE.eq(ObjectXrefEnsemblObjectType.Gene))
                    .fetch()
                    .stream()
                    .map(r -> r.get(XREF.DBPRIMARY_ACC))
                    .collect(Collectors.toList());

            final GeneAnnotation geneAnnotation = new GeneAnnotation(
                    variant, isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand);

            final Result<?> transcripts = mEnsemblDAO.queryTranscripts(geneId);

            for (final Record transcriptRecord : transcripts)
            {
                Transcript transcript = mEnsemblDAO.buildTranscript(
                        geneAnnotation, transcriptRecord, position, canonicalTranscriptId, geneStrand > 0);

                if (transcript != null)
                    geneAnnotation.addTranscript(transcript);
            }

            if (!geneAnnotation.transcripts().isEmpty())
                result.add(geneAnnotation);
        }

        return result;
    }

}

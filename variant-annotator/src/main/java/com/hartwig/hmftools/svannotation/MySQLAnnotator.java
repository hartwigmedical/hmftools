package com.hartwig.hmftools.svannotation;

import static java.lang.Math.abs;

import static org.ensembl.database.homo_sapiens_core.tables.CoordSystem.COORD_SYSTEM;
import static org.ensembl.database.homo_sapiens_core.tables.Exon.EXON;
import static org.ensembl.database.homo_sapiens_core.tables.ExonTranscript.EXON_TRANSCRIPT;
import static org.ensembl.database.homo_sapiens_core.tables.Gene.GENE;
import static org.ensembl.database.homo_sapiens_core.tables.Karyotype.KARYOTYPE;
import static org.ensembl.database.homo_sapiens_core.tables.ObjectXref.OBJECT_XREF;
import static org.ensembl.database.homo_sapiens_core.tables.SeqRegion.SEQ_REGION;
import static org.ensembl.database.homo_sapiens_core.tables.Transcript.TRANSCRIPT;
import static org.ensembl.database.homo_sapiens_core.tables.Translation.TRANSLATION;
import static org.ensembl.database.homo_sapiens_core.tables.Xref.XREF;
import static org.jooq.impl.DSL.boolAnd;
import static org.jooq.impl.DSL.decode;
import static org.jooq.impl.DSL.groupConcatDistinct;
import static org.jooq.impl.DSL.when;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.StructuralVariantAnnotation;
import com.hartwig.hmftools.common.variant.structural.annotation.Transcript;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.ensembl.database.homo_sapiens_core.enums.GeneStatus;
import org.ensembl.database.homo_sapiens_core.enums.ObjectXrefEnsemblObjectType;
import org.ensembl.database.homo_sapiens_core.tables.Exon;
import org.ensembl.database.homo_sapiens_core.tables.Xref;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.Condition;
import org.jooq.DSLContext;
import org.jooq.Query;
import org.jooq.Record;
import org.jooq.Result;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;
import org.jooq.types.UInteger;
import org.jooq.types.ULong;

public class MySQLAnnotator implements VariantAnnotator
{
    private static final String ENTREZ_IDS = "ENTREZ_IDS";
    private static final String KARYOTYPE_BAND = "KARYOTYPE_BAND";
    private static final String CODING_START = "CODING_START";
    private static final String CODING_END = "CODING_END";

    @NotNull
    private final DSLContext context;
    @NotNull
    private final UInteger coordSystemId;

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
        context = DSL.using(conn, SQLDialect.MYSQL);
        coordSystemId = findCoordSystemId();
    }

    public MySQLAnnotator(@NotNull final DSLContext dbConnection)
    {
        context = dbConnection;
        coordSystemId = findCoordSystemId();
    }

    @NotNull
    private UInteger findCoordSystemId()
    {
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
    public List<StructuralVariantAnnotation> annotateVariants(@NotNull List<EnrichedStructuralVariant> variants)
    {
        return variants.stream().map(this::annotateVariant).collect(Collectors.toList());
    }

    @NotNull
    private StructuralVariantAnnotation annotateVariant(@NotNull EnrichedStructuralVariant variant)
    {
        final StructuralVariantAnnotation annotation = new StructuralVariantAnnotation(variant);

        annotation.annotations().addAll(annotateBreakend(variant, true, variant.start().chromosome(), variant.start().position()));

        if (variant.end() != null)
        {
            annotation.annotations().addAll(annotateBreakend(variant, false, variant.end().chromosome(), variant.end().position()));
        }

        return annotation;
    }

    @NotNull
    private List<GeneAnnotation> annotateBreakend(@NotNull EnrichedStructuralVariant variant, final boolean isStart, @NotNull String chromosome, final long position)
    {
        final List<GeneAnnotation> result = Lists.newArrayList();

        final Result<?> genes = queryGenesOnChromosomeAndPosition(chromosome, position);

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

            final GeneAnnotation geneAnnotation = new GeneAnnotation(variant, isStart, geneName, geneStableId, geneStrand, synonyms, entrezIds, karyotypeBand);

            final Result<?> transcripts = queryTranscripts(geneId);

            for (final Record transcriptRecord : transcripts)
            {
                Transcript transcript = buildTranscript(geneAnnotation, transcriptRecord, position, canonicalTranscriptId, geneStrand > 0);
                // Transcript transcript = buildTranscript_old(geneAnnotation, transcriptRecord, position, canonicalTranscriptId, geneStrand > 0);

                if (transcript != null)
                    geneAnnotation.addTranscript(transcript);
            }

            if (!geneAnnotation.transcripts().isEmpty())
                result.add(geneAnnotation);
        }

        return result;
    }

    @NotNull
    private Result<?> queryGenesOnChromosomeAndPosition(@NotNull String chromosome, long position)
    {
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

    @NotNull
    private Result<?> queryTranscripts(@NotNull final UInteger geneId)
    {
        // CHSH: query created manually since the structured version below doens't produce the same result
        final String queryStr = "select t.transcript_id, t.stable_id,"
                + " if(t.seq_region_strand = -1, ce.seq_region_end - tl.seq_end + 1, cs.seq_region_start + tl.seq_start - 1) as CODING_START,"
                + " if(t.seq_region_strand = -1, cs.seq_region_end - tl.seq_start + 1, ce.seq_region_start + tl.seq_end - 1) as CODING_END"
                + " from transcript as t"
                + " left join translation tl on tl.transcript_id = t.transcript_id"
                + " left join exon cs on cs.exon_id = tl.start_exon_id"
                + " left join exon ce on ce.exon_id = tl.end_exon_id"
                + " where t.gene_id = '" + geneId + "'";

        return context.fetch(queryStr);

        /*
        final Exon EXON_START = EXON.as("cs");
        final Exon EXON_END = EXON.as("ce");

        return context.select(TRANSCRIPT.TRANSCRIPT_ID,
                TRANSCRIPT.STABLE_ID,
                when(TRANSCRIPT.SEQ_REGION_STRAND.eq((byte) -1), EXON_END.SEQ_REGION_END.minus(TRANSLATION.SEQ_END).plus(1)).otherwise(
                        EXON_START.SEQ_REGION_START).plus(TRANSLATION.SEQ_START).minus(1).as(CODING_START),
                when(TRANSCRIPT.SEQ_REGION_STRAND.eq((byte) -1), EXON_START.SEQ_REGION_END.minus(TRANSLATION.SEQ_START).plus(1)).otherwise(
                        EXON_END.SEQ_REGION_START).plus(TRANSLATION.SEQ_END).minus(1).as(CODING_END))
                .from(TRANSCRIPT)
                .leftJoin(TRANSLATION)
                .on(TRANSLATION.TRANSCRIPT_ID.eq(TRANSCRIPT.TRANSCRIPT_ID))
                .leftJoin(EXON_START)
                .on(EXON_START.EXON_ID.eq(TRANSLATION.START_EXON_ID))
                .leftJoin(EXON_END)
                .on(EXON_END.EXON_ID.eq(TRANSLATION.END_EXON_ID))
                .where(TRANSCRIPT.GENE_ID.eq(geneId))
                .fetch();
        */

    }

    private static String SPECIFIC_TRANS = "";
    // private static String SPECIFIC_TRANS = "ENST00000552327";

    private Transcript buildTranscript(
            @NotNull GeneAnnotation parent, @NotNull Record transcript, long position,
            @NotNull UInteger canonicalTranscriptId, boolean isForwardStrand)
    {
        // get all exons for this position and use them to work out:
        // exon up and downstream, and phasing
        // coding bases since start
        // total coding bases
        // total exons
        final UInteger transcriptId = transcript.get(TRANSCRIPT.TRANSCRIPT_ID);
        final boolean canonical = transcriptId.equals(canonicalTranscriptId);
        final String transcriptStableId = transcript.get(TRANSCRIPT.STABLE_ID);

        final Result<?> allExons = context.select(EXON_TRANSCRIPT.RANK, EXON.PHASE, EXON.END_PHASE, EXON.SEQ_REGION_START, EXON.SEQ_REGION_END)
                .from(EXON_TRANSCRIPT)
                .innerJoin(EXON)
                .on(EXON.EXON_ID.eq(EXON_TRANSCRIPT.EXON_ID))
                .where(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(transcriptId))
                .orderBy(EXON.SEQ_REGION_START.asc())
                .fetch();

        if(transcriptStableId.equals(SPECIFIC_TRANS))
        {
            LOGGER.debug("specific trans: {}", transcriptStableId);
        }

        int exonMax = allExons.size();

        ULong codingStartVal = (ULong) transcript.get(CODING_START);
        ULong codingEndVal = (ULong) transcript.get(CODING_END);

        boolean isCoding = codingStartVal != null && codingEndVal != null;
        long codingStart = codingStartVal != null ? codingStartVal.longValue() : 0;
        long codingEnd =  codingEndVal != null ? codingEndVal.longValue() : 0;

        // for the given position, determine how many coding bases occur prior to the position
        // in the direction of the transcript
        // strand direction will be corrected for afterwards

        boolean inCodingRegion = false;
        boolean codingRegionEnded = false;

        long codingBases = 0;
        long totalCodingBases = 0;
        int prevExonRank = -1;
        int prevExonPhase = 0;
        int prevExonEndPhase = 0;
        int nextExonRank = -1;
        int nextExonPhase = 0;
        int nextExonEndPhase = 0;

        for (final Record exon : allExons)
        {
            long exonStart = exon.get(EXON.SEQ_REGION_START).longValue();
            long exonEnd = exon.get(EXON.SEQ_REGION_END).longValue();

            if(position >= exonStart && position <= exonEnd)
            {
                // falls within an exon
                prevExonRank = exon.get(EXON_TRANSCRIPT.RANK);
                prevExonPhase = exon.get(EXON.PHASE);
                prevExonEndPhase = exon.get(EXON.END_PHASE);

                nextExonRank = exon.get(EXON_TRANSCRIPT.RANK);
                nextExonPhase = exon.get(EXON.PHASE);
                nextExonEndPhase = exon.get(EXON.END_PHASE);
            }
            else if(position > exonEnd)
            {
                // continue updating this until past the position
                prevExonRank = exon.get(EXON_TRANSCRIPT.RANK);
                prevExonPhase = exon.get(EXON.PHASE);
                prevExonEndPhase = exon.get(EXON.END_PHASE);
            }
            else if(position < exonStart && nextExonRank == -1)
            {
                // set at the first exon past this position
                nextExonRank = exon.get(EXON_TRANSCRIPT.RANK);
                nextExonPhase = exon.get(EXON.PHASE);
                nextExonEndPhase = exon.get(EXON.END_PHASE);
            }

            if(!isCoding)
                continue;

            if(!inCodingRegion)
            {
                if(exonEnd >= codingStart)
                {
                    // coding region begins in this exon
                    inCodingRegion = true;

                    totalCodingBases += exonEnd - codingStart + 1;

                    // check whether the position falls in this exon and if so before or after the coding start
                    if(position >= codingStart)
                    {
                        if(position < exonEnd)
                            codingBases += position - codingStart + 1;
                        else
                            codingBases += exonEnd - codingStart + 1;
                    }
                }
            }
            else if(!codingRegionEnded)
            {
                if(exonStart > codingEnd)
                {
                    codingRegionEnded = true;
                }
                else if(exonEnd > codingEnd)
                {
                    // coding region ends in this exon
                    codingRegionEnded = true;

                    totalCodingBases += codingEnd - exonStart + 1;

                    if(position >= exonStart)
                    {
                        if (position < codingEnd)
                            codingBases += position - exonStart + 1;
                        else
                            codingBases += codingEnd - exonStart + 1;
                    }
                }
                else
                {
                    // take all of the exon's bases
                    totalCodingBases += exonEnd - exonStart + 1;

                    if(position >= exonStart)
                    {
                        if (position < exonEnd)
                            codingBases += position - exonStart + 1;
                        else
                            codingBases += exonEnd - exonStart + 1;
                    }
                }
            }

            /*
            LOGGER.debug("transcript({}:{}) exon({}: {} - {}) coding({} - {}) position({}) coding(pos={} total={}) region(coding={} ended={})",
                    transcriptId, transcriptStableId, exon.get(EXON_TRANSCRIPT.RANK), exonStart, exonEnd,
                    codingStart, codingEnd, position, codingBases, totalCodingBases, inCodingRegion, codingRegionEnded);
            */

            if(codingBases < 0 || totalCodingBases < 0)
            {
                LOGGER.warn("transcript({}:{}) exon({}: {} - {}) coding({} - {}) position({}) coding(pos={} total={}) region(coding={} ended={})",
                        transcriptId, transcriptStableId, exon.get(EXON_TRANSCRIPT.RANK), exonStart, exonEnd,
                        codingStart, codingEnd, position, codingBases, totalCodingBases, inCodingRegion, codingRegionEnded);
                break;
            }
        }

        if(prevExonRank == -1)
        {
            if(!isForwardStrand)
            {
                // falls after the last exon on forward strand or before the first on reverse strand makes this position downstream
                LOGGER.debug("skipping transcript({}) position({}) after exon rank({} vs max={}) on reverse strand",
                        transcriptStableId, position, nextExonRank, exonMax);
                return null;
            }
            else
            {
                prevExonRank = 0;
                prevExonPhase = -1;
                prevExonEndPhase = -1;
            }
        }
        else if(nextExonRank == -1)
        {
            if(isForwardStrand)
            {
                // falls after the last exon on forward strand or before the first on reverse strand makes this position downstream
                LOGGER.debug("skipping transcript({}) position({}) after exon rank({} vs max={}) on forward strand",
                        transcriptStableId, position, prevExonRank, exonMax);
                return null;
            }
            else
            {
                nextExonRank = 0;
                nextExonPhase = -1;
                nextExonEndPhase = -1;
            }
        }

        if(nextExonRank < 0 || prevExonRank < 0 || abs(nextExonRank - prevExonRank) > 1)
        {
            LOGGER.warn("transcript({}) invalid exon ranks(prev={} next={}) forwardStrand({}) position({})",
                    transcriptStableId, prevExonRank, nextExonRank, isForwardStrand, position);
            return null;
        }

        if(isForwardStrand)
        {
            return new Transcript(parent,
                    transcriptStableId,
                    prevExonRank, prevExonEndPhase,
                    nextExonRank, nextExonPhase,
                    codingBases, totalCodingBases,
                    exonMax, canonical,
                    codingStartVal != null ? codingStart : null,
                    codingEndVal != null ? codingEnd : null);
        }
        else
        {
            return new Transcript(parent,
                    transcriptStableId,
                    nextExonRank, nextExonEndPhase, // note the switch
                    prevExonRank, prevExonPhase,
                    totalCodingBases - codingBases, totalCodingBases,
                    exonMax, canonical,
                    codingStartVal != null ? codingStart : null,
                    codingEndVal != null ? codingEnd : null);
        }
    }

    @Nullable
    private Transcript buildTranscript_old(
            @NotNull GeneAnnotation parent, @NotNull Record transcript, long position,
            @NotNull UInteger canonicalTranscriptId, boolean isForwardStrand)
    {
        final UInteger transcriptId = transcript.get(TRANSCRIPT.TRANSCRIPT_ID);
        final boolean canonical = transcriptId.equals(canonicalTranscriptId);
        final String transcriptStableId = transcript.get(TRANSCRIPT.STABLE_ID);

        final Record exonLeft = context.select(EXON_TRANSCRIPT.RANK, EXON.PHASE, EXON.END_PHASE, EXON.SEQ_REGION_START, EXON.SEQ_REGION_END)
                .from(EXON_TRANSCRIPT)
                .innerJoin(EXON)
                .on(EXON.EXON_ID.eq(EXON_TRANSCRIPT.EXON_ID))
                .where(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(transcriptId))
                .and(EXON.SEQ_REGION_START.le(UInteger.valueOf(position)))
                .orderBy(EXON.SEQ_REGION_START.desc())
                .limit(1)
                .fetchOne();

        final Record exonRight = context.select(EXON_TRANSCRIPT.RANK, EXON.PHASE, EXON.END_PHASE, EXON.SEQ_REGION_START, EXON.SEQ_REGION_END)
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

        if (isForwardStrand)
        {
            exonUpstream = exonLeft == null ? 0 : exonLeft.get(EXON_TRANSCRIPT.RANK);
            exonUpstreamPhase = exonLeft == null ? -1 : exonLeft.get(EXON.END_PHASE);

            exonDownstream = exonRight == null ? 0 : exonRight.get(EXON_TRANSCRIPT.RANK);
            exonDownstreamPhase = exonRight == null ? -1 : exonRight.get(EXON.PHASE);
        }
        else
        {
            exonDownstream = exonLeft == null ? 0 : exonLeft.get(EXON_TRANSCRIPT.RANK);
            exonDownstreamPhase = exonLeft == null ? -1 : exonLeft.get(EXON.PHASE);

            exonUpstream = exonRight == null ? 0 : exonRight.get(EXON_TRANSCRIPT.RANK);
            exonUpstreamPhase = exonRight == null ? -1 : exonRight.get(EXON.END_PHASE);
        }

        if (exonUpstream > 0 && exonDownstream == 0)
        {
            // past the last exon
            return null;
        }
        else
        {
            UInteger codingStart = (UInteger) transcript.get(CODING_START);
            UInteger codingEnd = (UInteger) transcript.get(CODING_END);

            long exonStart = 0;
            long exonEnd = 0;

            if(exonUpstream == exonDownstream && exonRight != null)
            {
                exonStart = ((UInteger)exonRight.get(EXON.SEQ_REGION_START)).longValue();
                exonEnd = ((UInteger)exonRight.get(EXON.SEQ_REGION_END)).longValue();
            }

            return new Transcript(parent,
                    transcriptStableId,
                    exonUpstream,
                    exonUpstreamPhase,
                    exonDownstream,
                    exonDownstreamPhase,
                    exonStart,
                    exonEnd,
                    exonMax,
                    canonical,
                    codingStart != null ? codingStart.longValue() : null,
                    codingEnd != null ? codingEnd.longValue() : null);
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

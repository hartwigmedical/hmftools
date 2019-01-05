package com.hartwig.hmftools.svannotation;

import static java.lang.Integer.max;
import static java.lang.Math.abs;
import static java.lang.Math.min;

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
import com.hartwig.hmftools.svannotation.analysis.RnaFusionData;

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
import org.jooq.types.ULong;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EnsemblDAO
{
    private static final String ENTREZ_IDS = "ENTREZ_IDS";
    private static final String KARYOTYPE_BAND = "KARYOTYPE_BAND";
    private static final String CODING_START = "CODING_START";
    private static final String CODING_END = "CODING_END";

    @NotNull
    private final DSLContext context;
    @NotNull
    private final UInteger coordSystemId;

    private static final Logger LOGGER = LogManager.getLogger(EnsemblDAO.class);

    public EnsemblDAO(@NotNull final DSLContext dbConnection)
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

    @NotNull
    public Result<?> queryGenesOnChromosomeAndPosition(@NotNull String chromosome, long position, int preGeneDistance)
    {
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
                        decode().when(GENE.SEQ_REGION_START.ge(UInteger.valueOf(preGeneDistance)),
                                GENE.SEQ_REGION_START.sub(preGeneDistance)).otherwise(GENE.SEQ_REGION_START))
                        .otherwise(GENE.SEQ_REGION_START)
                        .le(UInteger.valueOf(position)))
                .and(decode().when(GENE.SEQ_REGION_STRAND.lt(zero), GENE.SEQ_REGION_END.add(preGeneDistance))
                        .otherwise(GENE.SEQ_REGION_END)
                        .ge(UInteger.valueOf(position)))
                .and(geneStartInKaryotypeBand().or(geneEndInKaryotypeBand()))
                .groupBy(GENE.GENE_ID)
                .fetch();
    }

    @NotNull
    public Result<?> queryTranscripts(@NotNull final UInteger geneId)
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

    // private static String SPEC_TRANSCRIPT = "ENST00000453137";
    private static String SPEC_TRANSCRIPT = "";

    public Transcript buildTranscript(
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

        /*
        if(transcriptStableId.equals(SPEC_TRANSCRIPT))
        {
            LOGGER.debug("specific transcript({})", transcriptStableId);
        }
        */

        final Result<?> allExons = context.select(EXON_TRANSCRIPT.RANK, EXON.PHASE, EXON.END_PHASE, EXON.SEQ_REGION_START, EXON.SEQ_REGION_END)
                .from(EXON_TRANSCRIPT)
                .innerJoin(EXON)
                .on(EXON.EXON_ID.eq(EXON_TRANSCRIPT.EXON_ID))
                .where(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(transcriptId))
                .orderBy(EXON.SEQ_REGION_START.asc())
                .fetch();

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
        long transcriptStart = 0;
        long transcriptEnd = 0;

        // previous here will be the earlier exon, ordered by increasing position (ie regardless of strand direction)
        int prevExonRank = -1;
        int prevExonPhase = 0;
        int prevExonEndPhase = 0;

        // similarly the next exon will be exon immediately after the position for exons which increase with positino
        int nextExonRank = -1;
        int nextExonPhase = 0;
        int nextExonEndPhase = 0;

        for (int index = 0; index < allExons.size(); ++index)
        {
            final Record exon = allExons.get(index);

            long exonStart = exon.get(EXON.SEQ_REGION_START).longValue();
            long exonEnd = exon.get(EXON.SEQ_REGION_END).longValue();

            if(index == 0)
                transcriptStart = exonStart;

            if(index == allExons.size() - 1)
                transcriptEnd = exonEnd;

            if(position >= exonStart && position <= exonEnd)
            {
                // falls within an exon
                prevExonRank = nextExonRank = exon.get(EXON_TRANSCRIPT.RANK);

                if(isForwardStrand)
                {
                    prevExonEndPhase = exon.get(EXON.PHASE);
                    nextExonPhase = exon.get(EXON.END_PHASE);

                    // won't be used
                    prevExonPhase = prevExonEndPhase;
                    nextExonEndPhase = nextExonPhase;
                }
                else
                {
                    prevExonPhase = exon.get(EXON.END_PHASE);
                    nextExonEndPhase = exon.get(EXON.PHASE);

                    prevExonEndPhase = prevExonPhase;
                    nextExonPhase = nextExonEndPhase;
                }
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
                //                LOGGER.debug("skipping transcript({}) position({}) after exon rank({} vs max={}) on reverse strand",
                //                        transcriptStableId, position, nextExonRank, exonMax);
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
                //                LOGGER.debug("skipping transcript({}) position({}) after exon rank({} vs max={}) on forward strand",
                //                        transcriptStableId, position, prevExonRank, exonMax);
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
                    exonMax, canonical, transcriptStart, transcriptEnd,
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
                    exonMax, canonical, transcriptStart, transcriptEnd,
                    codingStartVal != null ? codingStart : null,
                    codingEndVal != null ? codingEnd : null);
        }
    }

    public static int EXON_DATA_TRANSCRIPT_ID = 0;
    public static int EXON_DATA_EXON_RANK = 1;

    public String getExonDetailsForPosition(final String chromosome, long posStart, long posEnd)
    {
        /*
        select t.transcript_id, t.stable_id, e.seq_region_start, e.seq_region_end, e.seq_region_end - e.seq_region_start as exonLength, et.rank, x.display_label
        from transcript as t, exon as e, exon_transcript as et, gene as g, xref as x, seq_region as sr
        where t.transcript_id = et.transcript_id and e.exon_id = et.exon_id and g.gene_id = t.gene_id
        and g.gene_id = t.gene_id and sr.seq_region_id = g.seq_region_id
        and sr.name = 11
        and g.display_xref_id = x.xref_id
        and 4127276 >= e.seq_region_start -10 && 4127276 <= e.seq_region_end + 10;
         */

        // first just look for matching exons prior to any concern for their transcripts
        final Result<?> exonDataList = context.select(EXON.EXON_ID, EXON.SEQ_REGION_START, EXON.SEQ_REGION_END)
                .from(EXON, SEQ_REGION)
                .where(EXON.SEQ_REGION_ID.eq(SEQ_REGION.SEQ_REGION_ID))
                .and(EXON.SEQ_REGION_START.greaterOrEqual(UInteger.valueOf(posStart - 1)))
                .and(EXON.SEQ_REGION_START.lessOrEqual(UInteger.valueOf(posStart + 1)))
                .and(EXON.SEQ_REGION_END.greaterOrEqual(UInteger.valueOf(posEnd - 1)))
                .and(EXON.SEQ_REGION_END.lessOrEqual(UInteger.valueOf(posEnd + 1)))
                .fetch();

        if(exonDataList.isEmpty())
            return "";

        String exonDataStr = "";

        for (final Record exonData : exonDataList)
        {
            final UInteger exonId = exonData.get(EXON.EXON_ID);
            long exonStart = exonData.get(EXON.SEQ_REGION_START).longValue();
            long exonEnd = exonData.get(EXON.SEQ_REGION_END).longValue();

            // now retrieve transcript info
            final Result<?> transcriptDataList = context.select(TRANSCRIPT.TRANSCRIPT_ID, TRANSCRIPT.STABLE_ID,
                    EXON_TRANSCRIPT.RANK, GENE.CANONICAL_TRANSCRIPT_ID)
                    .from(GENE, EXON_TRANSCRIPT, TRANSCRIPT)
                    .where(TRANSCRIPT.GENE_ID.eq(GENE.GENE_ID))
                    .and(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(TRANSCRIPT.TRANSCRIPT_ID))
                    .and(EXON_TRANSCRIPT.EXON_ID.eq(exonId))
                    .fetch();

            for(final Record transcriptData : transcriptDataList)
            {
                int exonRank = transcriptData.get(EXON_TRANSCRIPT.RANK);
                final String transcriptStableId = transcriptData.get(TRANSCRIPT.STABLE_ID);
                final UInteger transcriptId = transcriptData.get(TRANSCRIPT.TRANSCRIPT_ID);
                final UInteger canonicalTranscriptId = transcriptData.get(GENE.CANONICAL_TRANSCRIPT_ID);

                if(exonDataStr.isEmpty() || canonicalTranscriptId == transcriptId)
                {
                    exonDataStr = String.format("%s;%d;%d", transcriptStableId, exonRank, exonEnd - exonStart);

                    if(canonicalTranscriptId == transcriptId)
                        break;
                }
            }
        }

        return exonDataStr;
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

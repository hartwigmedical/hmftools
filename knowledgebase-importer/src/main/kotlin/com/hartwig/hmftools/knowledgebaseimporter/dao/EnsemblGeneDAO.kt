package com.hartwig.hmftools.knowledgebaseimporter.dao

import org.ensembl.database.homo_sapiens_core.Tables.*
import org.jooq.*
import org.jooq.impl.DSL
import org.jooq.types.UInteger
import java.sql.DriverManager
import java.util.stream.Stream
import kotlin.streams.toList

typealias ExonRecord = Record5<UInteger, UInteger, Byte, Byte, Byte>
typealias TranslationRecord = Record2<Int, Int>

class EnsemblGeneDAO(url: String, user: String, password: String) {
    companion object {
        private fun baseExonQuery(context: DSLContext): SelectOnConditionStep<ExonRecord> {
            return context.select(EXON.SEQ_REGION_START, EXON.SEQ_REGION_END, EXON.SEQ_REGION_STRAND, EXON.PHASE, EXON.END_PHASE)
                    .from(EXON)
                    .join(EXON_TRANSCRIPT).on(EXON.EXON_ID.eq(EXON_TRANSCRIPT.EXON_ID))
        }

        private fun baseTranslationQuery(context: DSLContext): SelectJoinStep<TranslationRecord> {
            return context.select(TRANSLATION.SEQ_START, TRANSLATION.SEQ_END).from(TRANSLATION)
        }

        private val createExonsLambda: (exonRecords: Stream<ExonRecord>) -> List<Exon> = {
            it.toList().map {
                val strand = it[EXON.SEQ_REGION_STRAND].toInt()
                val start = it[EXON.SEQ_REGION_START].toLong()
                val end = it[EXON.SEQ_REGION_END].toLong()
                if (strand > 0) {
                    Exon(start, end, it[EXON.PHASE].toInt(), it[EXON.END_PHASE].toInt())
                } else {
                    Exon(end * -1, start * -1, it[EXON.PHASE].toInt(), it[EXON.END_PHASE].toInt())
                }
            }
        }

        private val createGeneStartEndLambda: (translationRecords: Stream<TranslationRecord>) -> Pair<Long, Long> = {
            val translationList = it.toList().map { Pair(it[TRANSLATION.SEQ_START].toLong(), it[TRANSLATION.SEQ_END].toLong()) }
            assert(translationList.size == 1)
            translationList.first()
        }
    }

    private val context = DSL.using(DriverManager.getConnection(url, user, password), SQLDialect.MYSQL)

    fun canonicalGeneModel(geneName: String): Gene {
        val exons = canonicalExons(geneName)
        val (seqStart, seqEnd) = canonicalSequenceStartEnd(geneName)
        return Gene(exons, seqStart, seqEnd)
    }

    fun transcriptGeneModel(transcript: String): Gene {
        val exons = transcriptExons(transcript)
        val (seqStart, seqEnd) = transcriptSequenceStartEnd(transcript)
        return Gene(exons, seqStart, seqEnd)
    }

    private fun canonicalExons(geneName: String): List<Exon> {
        return baseExonQuery(context).run {
            this.join(GENE).on(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(GENE.CANONICAL_TRANSCRIPT_ID))
                    .join(XREF).on(XREF.XREF_ID.eq(GENE.DISPLAY_XREF_ID))
                    .where(XREF.DISPLAY_LABEL.eq(geneName))
                    .stream()
        }.run(createExonsLambda)
    }

    private fun canonicalSequenceStartEnd(geneName: String): Pair<Long, Long> {
        return baseTranslationQuery(context).run {
            this.join(GENE).on(TRANSLATION.TRANSCRIPT_ID.eq(GENE.CANONICAL_TRANSCRIPT_ID))
                    .join(XREF).on(XREF.XREF_ID.eq(GENE.DISPLAY_XREF_ID))
                    .where(XREF.DISPLAY_LABEL.eq(geneName))
                    .stream()
        }.run(createGeneStartEndLambda)
    }


    private fun transcriptExons(transcript: String): List<Exon> {
        return baseExonQuery(context).run {
            this.join(TRANSCRIPT).on(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(TRANSCRIPT.TRANSCRIPT_ID))
                    .where(TRANSCRIPT.STABLE_ID.eq(transcript))
                    .stream()
        }.run(createExonsLambda)
    }

    private fun transcriptSequenceStartEnd(transcript: String): Pair<Long, Long> {
        return baseTranslationQuery(context).run {
            this.join(TRANSCRIPT).on(TRANSLATION.TRANSCRIPT_ID.eq(TRANSCRIPT.TRANSCRIPT_ID))
                    .where(TRANSCRIPT.STABLE_ID.eq(transcript))
                    .stream()
        }.run(createGeneStartEndLambda)
    }
}

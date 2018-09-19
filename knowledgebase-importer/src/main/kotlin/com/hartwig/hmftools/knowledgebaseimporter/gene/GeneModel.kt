package com.hartwig.hmftools.knowledgebaseimporter.gene

import com.hartwig.hmftools.common.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.common.region.HmfTranscriptRegion
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GenericMutation
import org.apache.logging.log4j.LogManager
import org.jooq.Record2
import org.jooq.Record6
import org.jooq.SQLDialect
import org.jooq.impl.DSL
import org.jooq.types.UInteger
import java.sql.DriverManager

typealias ExonRecord = Record6<String, UInteger, UInteger, Byte, Byte, Byte>
typealias TranslationRecord = Record2<Int, Int>

class GeneModel(ensemblDbUrl: String, hmfpatientsDbUrl: String, user: String, password: String) {
    companion object {
        private val hmfGenePanel by lazy { HmfGenePanelSupplier.allGenesMap() }
        private val logger = LogManager.getLogger("GeneModel")

//        private fun baseExonQuery(context: DSLContext): SelectOnConditionStep<ExonRecord> {
//            return context.select(SEQ_REGION.NAME, EXON.SEQ_REGION_START, EXON.SEQ_REGION_END, EXON.SEQ_REGION_STRAND, EXON.PHASE,
//                    EXON.END_PHASE)
//                    .from(EXON)
//                    .join(EXON_TRANSCRIPT).on(EXON.EXON_ID.eq(EXON_TRANSCRIPT.EXON_ID))
//                    .join(SEQ_REGION).on(EXON.SEQ_REGION_ID.eq(SEQ_REGION.SEQ_REGION_ID))
//        }
//
//        private fun startExonQuery(context: DSLContext): SelectOnConditionStep<ExonRecord> {
//            return baseExonQuery(context).join(TRANSLATION).on(TRANSLATION.START_EXON_ID.eq(EXON.EXON_ID))
//        }
//
//        private fun endExonQuery(context: DSLContext): SelectOnConditionStep<ExonRecord> {
//            return baseExonQuery(context).join(TRANSLATION).on(TRANSLATION.END_EXON_ID.eq(EXON.EXON_ID))
//        }
//
//        private fun baseTranslationQuery(context: DSLContext): SelectJoinStep<TranslationRecord> {
//            return context.select(TRANSLATION.SEQ_START, TRANSLATION.SEQ_END).from(TRANSLATION)
//        }
//
//        private fun canonicalTranscriptQuery(hmfContext: DSLContext, geneName: String): List<String> {
//            return hmfContext.select(CANONICALTRANSCRIPT.TRANSCRIPTID).from(CANONICALTRANSCRIPT)
//                    .where(CANONICALTRANSCRIPT.GENE.eq(geneName)).stream().map { it.get(CANONICALTRANSCRIPT.TRANSCRIPTID) }.toList()
//        }
//
//        private val createExonsLambda: (exonRecords: Stream<ExonRecord>) -> List<Exon> = {
//            it.toList().map {
//                val chromosome = it[SEQ_REGION.NAME]
//                val strand = it[EXON.SEQ_REGION_STRAND].toInt()
//                val start = it[EXON.SEQ_REGION_START].toLong()
//                val end = it[EXON.SEQ_REGION_END].toLong()
//                if (strand > 0) {
//                    Exon(chromosome, start, end, it[EXON.PHASE].toInt(), it[EXON.END_PHASE].toInt())
//                } else {
//                    Exon(chromosome, end * -1, start * -1, it[EXON.PHASE].toInt(), it[EXON.END_PHASE].toInt())
//                }
//            }
//        }
//
//        private fun createGeneStartEndLambda(translationRecords: Stream<TranslationRecord>, gene: String): Pair<Long, Long>? {
//            val translationList = translationRecords.toList().map {
//                Pair(it[TRANSLATION.SEQ_START].toLong(), it[TRANSLATION.SEQ_END].toLong())
//            }
//            if (translationList.size > 1) logger.warn("Expected max 1 translation for gene $gene but found ${translationList.size}")
//            return translationList.firstOrNull()
//        }
    }

    private val hmfContext = DSL.using(DriverManager.getConnection(hmfpatientsDbUrl, user, password), SQLDialect.MYSQL)
    private val context = DSL.using(DriverManager.getConnection(ensemblDbUrl, user, password), SQLDialect.MYSQL)

    fun hmfTranscriptForGene(gene: String): String? {
        // TODO (KODU): Map MLL2 gene to KMT2B
        val transcriptRegion = hmfGenePanel[gene]
        if (transcriptRegion == null) {
            logger.warn("Gene $gene not found in HMF gene panel!")
            return null
        }
        return transcriptRegion.transcriptID()
    }

    fun hmfTranscriptRegionForGenericMutation(mutation: GenericMutation): HmfTranscriptRegion? {
        // TODO (KODU): Map MLL2 gene to KMT2B
        val transcriptRegion = hmfGenePanel[mutation.gene]
        if (transcriptRegion == null) {
            logger.warn("Gene ${mutation.gene} not found in HMF gene panel!")
        } else if (mutation.transcript != null) {
            if (transcriptRegion.transcriptID() != mutation.transcript) {
                logger.warn("Non-canonical transcript ${mutation.transcript} requested for gene ${mutation.gene}")
                return null
            }
        }

        return transcriptRegion
    }

//    fun hmfCanonicalTranscript(geneName: String): List<String> {
//        val transcripts = canonicalTranscriptQuery(hmfContext, geneName)
//        if (transcripts.size != 1) logger.warn("Expected single canonical transcript for gene $geneName, but found ${transcripts.size}")
//        return transcripts
//    }
//
//
//    fun transcriptGeneModel(transcript: String): Gene? {
//        val exons = transcriptExons(transcript)
//        return transcriptSequenceStartEnd(transcript)?.run {
//            val startEnd = transcriptStartEnd(transcript).toList()
//            Gene(exons, startExon(transcript), endExon(transcript), first, second, startEnd.min()!!, startEnd.max()!!, transcript)
//        }
//    }
//
//    private fun canonicalExons(geneName: String): List<Exon> {
//        return baseExonQuery(context).join(GENE).on(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(GENE.CANONICAL_TRANSCRIPT_ID))
//                .join(XREF).on(XREF.XREF_ID.eq(GENE.DISPLAY_XREF_ID))
//                .where(XREF.DISPLAY_LABEL.eq(geneName))
//                .stream()
//                .run(createExonsLambda)
//    }
//
//    private fun canonicalSequenceStartEnd(geneName: String): Pair<Long, Long>? {
//        return baseTranslationQuery(context).join(GENE).on(TRANSLATION.TRANSCRIPT_ID.eq(GENE.CANONICAL_TRANSCRIPT_ID))
//                .join(XREF).on(XREF.XREF_ID.eq(GENE.DISPLAY_XREF_ID))
//                .where(XREF.DISPLAY_LABEL.eq(geneName))
//                .stream()
//                .run { createGeneStartEndLambda(this, geneName) }
//    }
//
//    private fun canonicalStartExon(geneName: String): Exon {
//        return startExonQuery(context).join(GENE).on(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(GENE.CANONICAL_TRANSCRIPT_ID))
//                .join(XREF).on(XREF.XREF_ID.eq(GENE.DISPLAY_XREF_ID))
//                .where(XREF.DISPLAY_LABEL.eq(geneName))
//                .stream()
//                .run(createExonsLambda)
//                .first()
//    }
//
//    private fun canonicalEndExon(geneName: String): Exon {
//        return endExonQuery(context).join(GENE).on(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(GENE.CANONICAL_TRANSCRIPT_ID))
//                .join(XREF).on(XREF.XREF_ID.eq(GENE.DISPLAY_XREF_ID))
//                .where(XREF.DISPLAY_LABEL.eq(geneName))
//                .stream()
//                .run(createExonsLambda)
//                .first()
//    }
//
//    private fun geneStartEnd(geneName: String): Pair<Long, Long> {
//        return context.select(GENE.SEQ_REGION_START, GENE.SEQ_REGION_END).from(GENE).join(XREF).on(XREF.XREF_ID.eq(GENE.DISPLAY_XREF_ID))
//                .where(XREF.DISPLAY_LABEL.eq(geneName))
//                .stream().toList()
//                .map { Pair(it[GENE.SEQ_REGION_START].toLong(), it[GENE.SEQ_REGION_END].toLong()) }.first()
//    }
//
//    private fun transcriptExons(transcript: String): List<Exon> {
//        return baseExonQuery(context).join(TRANSCRIPT).on(EXON_TRANSCRIPT.TRANSCRIPT_ID.eq(TRANSCRIPT.TRANSCRIPT_ID))
//                .where(TRANSCRIPT.STABLE_ID.eq(transcript))
//                .stream()
//                .run(createExonsLambda)
//    }

//    private fun transcriptSequenceStartEnd(transcript: String): Pair<Long, Long>? {
//        return baseTranslationQuery(context).join(TRANSCRIPT).on(TRANSLATION.TRANSCRIPT_ID.eq(TRANSCRIPT.TRANSCRIPT_ID))
//                .where(TRANSCRIPT.STABLE_ID.eq(transcript))
//                .stream()
//                .run { createGeneStartEndLambda(this, transcript) }
//    }
//
//    private fun startExon(transcript: String): Exon {
//        return startExonQuery(context)
//                .join(TRANSCRIPT).on(TRANSLATION.TRANSCRIPT_ID.eq(TRANSCRIPT.TRANSCRIPT_ID)).where(TRANSCRIPT.STABLE_ID.eq(transcript))
//                .stream()
//                .run(createExonsLambda)
//                .first()
//    }
//
//    private fun endExon(transcript: String): Exon {
//        return endExonQuery(context).join(TRANSCRIPT).on(TRANSLATION.TRANSCRIPT_ID.eq(TRANSCRIPT.TRANSCRIPT_ID))
//                .where(TRANSCRIPT.STABLE_ID.eq(transcript))
//                .stream()
//                .run(createExonsLambda)
//                .first()
//    }
//
//    private fun transcriptStartEnd(transcript: String): Pair<Long, Long> {
//        return context.select(TRANSCRIPT.SEQ_REGION_START, TRANSCRIPT.SEQ_REGION_END).from(TRANSCRIPT)
//                .where(TRANSCRIPT.STABLE_ID.eq(transcript))
//                .stream().toList()
//                .map { Pair(it[TRANSCRIPT.SEQ_REGION_START].toLong(), it[TRANSCRIPT.SEQ_REGION_END].toLong()) }.first()
//    }
}

package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.common.genome.bed.NamedBed
import com.hartwig.hmftools.common.genome.position.GenomePosition
import com.hartwig.hmftools.common.genome.position.GenomePositions
import com.hartwig.hmftools.common.genome.region.CodingRegions
import com.hartwig.hmftools.common.genome.region.GenomeRegions
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.variant.VariantContextDecorator
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideFragmentFactory
import com.hartwig.hmftools.lilac.sam.Indel
import com.hartwig.hmftools.lilac.sam.SAMCodingRecord
import com.hartwig.hmftools.lilac.sam.SAMSlicer
import htsjdk.samtools.SAMRecord
import htsjdk.samtools.SamReaderFactory
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.function.Consumer
import kotlin.math.abs

class SAMRecordReader(private val bamFile: String, private val refGenome: String, private val transcripts: List<HmfTranscriptRegion>, private val factory: NucleotideFragmentFactory) {

    companion object {
        const val MIN_MAPPING_QUALITY = 1
        const val MAX_DISTANCE = 1000
        val logger = LogManager.getLogger(this::class.java)
        val indelPon = SAMRecordReader::class.java.getResource("/pon/indels.txt")
                .readText()
                .split("\n")
                .map { Indel.fromString(it) }
                .toSet()

        val STOP_LOSS_ON_C = Indel("6", 31237115, "CN", "C")
    }

    private val codingRegions = transcripts.map { GenomeRegions.create(it.chromosome(), it.codingStart() - MAX_DISTANCE, it.codingEnd() + MAX_DISTANCE) }
    private val stopLossOnC = mutableMapOf<Indel, Int>()
    private val unmatchedIndels = mutableMapOf<Indel, Int>()
    private val unmatchedPONIndels = mutableMapOf<Indel, Int>()
    private var alignmentFiltered = 0

    fun alignmentFiltered(): Int {
        return alignmentFiltered
    }

    fun stopLossOnCIndels(): Int {
        return stopLossOnC[STOP_LOSS_ON_C] ?: 0
    }

    fun unmatchedIndels(minCount: Int): Map<Indel, Int> {
        return unmatchedIndels.filter { it.value >= minCount }
    }

    fun unmatchedPonIndels(minCount: Int): Map<Indel, Int> {
        return unmatchedPONIndels.filter { it.value >= minCount }
    }

    fun readFromBam(): List<NucleotideFragment> {
        return transcripts.flatMap { readFromBam(it, bamFile) }
    }

    fun readFromBam(variant: VariantContextDecorator): List<NucleotideFragment> {
        val variantPosition = GenomePositions.create(variant.chromosome(), variant.position())

        for (transcript in transcripts) {
            val reverseStrand = transcript.strand() == Strand.REVERSE
            val codingRegions = if (reverseStrand) codingRegions(transcript).reversed() else codingRegions(transcript)
            var hlaCodingRegionOffset = 0
            for (codingRegion in codingRegions) {
                if (codingRegion.contains(variantPosition) || abs(codingRegion.start() - variant.position()) <= 5 || abs(codingRegion.end() - variant.position()) <= 5) {
                    val codingRecords = query(reverseStrand, variantPosition, codingRegion, bamFile)
                            .filter { recordContainsVariant(variant, it) }

                    val distinctCodingRecords = codingRecords.distinct()

                    val nucleotideFragments = distinctCodingRecords
                            .mapNotNull { factory.createAlignmentFragments(it, codingRegion) }

                    val mateFragments = queryMateFragments(transcript, distinctCodingRecords)

                    return (nucleotideFragments + mateFragments)
                            .groupBy { it.id }
                            .map { it.value.reduce { x, y -> NucleotideFragment.merge(x, y) } }

                }
                hlaCodingRegionOffset += codingRegion.bases().toInt()
            }
        }

        return listOf()
    }

    private fun recordContainsVariant(variant: VariantContextDecorator, record: SAMCodingRecord): Boolean {
        if (variant.alt().length != variant.ref().length) {
            val expectedIndel = Indel(variant.chromosome(), variant.position().toInt(), variant.ref(), variant.alt())
            return record.indels.any { it.match(expectedIndel) }
        }

        for (i in variant.alt().indices) {
            val position = variant.position().toInt() + i
            val expectedBase = variant.alt()[i]
            val readIndex = record.record.getReadPositionAtReferencePosition(position) - 1
            if (readIndex < 0) {
                return false
            }

            if (record.record.readBases[readIndex].toChar() != expectedBase) {
                return false
            }

        }

        return true
    }


    private fun readFromBam(transcript: HmfTranscriptRegion, bamFile: String): List<NucleotideFragment> {
        logger.info("    querying ${transcript.gene()} (${transcript.chromosome()}:${transcript.codingStart()}-${transcript.codingEnd()})")

        val reverseStrand = transcript.strand() == Strand.REVERSE
        val codingRegions = if (reverseStrand) codingRegions(transcript).reversed() else codingRegions(transcript)

        val realignedRegions = mutableListOf<NucleotideFragment>()
        for (codingRegion in codingRegions) {
            realignedRegions.addAll(realign(codingRegion, reverseStrand, bamFile))
        }

        val groupedReads = realignedRegions.groupBy { it.id }
        val mergedReads = groupedReads.map { it.value.reduce { x, y -> NucleotideFragment.merge(x, y) } }

        return mergedReads
    }

    private fun codingRegions(transcript: HmfTranscriptRegion): List<NamedBed> {
        return CodingRegions.codingRegions(transcript)
    }

    private fun samReaderFactory(): SamReaderFactory {
        val default = SamReaderFactory.makeDefault()
        return if (refGenome.isNotEmpty()) {
            default.referenceSequence(File(refGenome))
        } else {
            default
        }
    }

    private fun queryMateFragments(transcript: HmfTranscriptRegion, codingRecords: List<SAMCodingRecord>): List<NucleotideFragment> {
        val slicer = SAMSlicer(MIN_MAPPING_QUALITY)
        val samRecords = codingRecords.map { it.record }.distinct()
        val mates = samReaderFactory().open(File(bamFile)).use { reader -> slicer.queryMates(reader, samRecords) }
        val result = mutableListOf<NucleotideFragment>()

        val reverseStrand = transcript.strand() == Strand.REVERSE
        val codingRegions = if (reverseStrand) codingRegions(transcript).reversed() else codingRegions(transcript)

        for (codingRegion in codingRegions) {
            mates
                    .filter { it.alignmentStart <= codingRegion.end() && it.alignmentEnd >= codingRegion.start() }
                    .map { SAMCodingRecord.create(reverseStrand, codingRegion, it) }
                    .mapNotNull { factory.createAlignmentFragments(it, codingRegion) }
                    .forEach { result.add(it) }
        }

        return result
    }


    private fun query(reverseStrand: Boolean, variantRegion: GenomePosition, nearestCodingRegion: NamedBed, bamFileName: String): List<SAMCodingRecord> {
        val slicer = SAMSlicer(MIN_MAPPING_QUALITY)
        val result = mutableListOf<SAMCodingRecord>()
        samReaderFactory().open(File(bamFileName)).use { samReader ->
            val consumer = Consumer<SAMRecord> { samRecord ->
                if (samRecord.bothEndsInRangeOfCodingTranscripts()) {
                    result.add(SAMCodingRecord.create(reverseStrand, nearestCodingRegion, samRecord))

                } else {
                    alignmentFiltered++
                }
            }

            slicer.slice(variantRegion.chromosome(), variantRegion.position().toInt(), variantRegion.position().toInt(), samReader, consumer)
        }
        return result
    }

    private fun query(reverseStrand: Boolean, codingRegion: NamedBed, bamFileName: String): List<SAMCodingRecord> {
        val slicer = SAMSlicer(MIN_MAPPING_QUALITY)
        val result = mutableListOf<SAMCodingRecord>()
        samReaderFactory().open(File(bamFileName)).use { samReader ->
            val consumer = Consumer<SAMRecord> { samRecord ->
                if (samRecord.bothEndsInRangeOfCodingTranscripts()) {
                    result.add(SAMCodingRecord.create(reverseStrand, codingRegion, samRecord))

                } else {
                    alignmentFiltered++
                }
            }

            slicer.slice(codingRegion, samReader, consumer)
        }
        return result
    }

    private fun realign(codingRegion: NamedBed, reverseStrand: Boolean, bamFileName: String): List<NucleotideFragment> {
        val result = mutableListOf<NucleotideFragment>()
        for (codingRecord in query(reverseStrand, codingRegion, bamFileName)) {
            if (codingRecord.indels.contains(STOP_LOSS_ON_C)) {
                stopLossOnC.compute(STOP_LOSS_ON_C) { _, u -> (u ?: 0) + 1 }
            }

            val fragment = factory.createFragment(codingRecord, codingRegion)
            if (fragment != null) {
                result.add(fragment)
            } else {
                for (indel in codingRecord.indels) {
                    if (indel in indelPon) {
                        unmatchedPONIndels.compute(indel) { _, u -> (u ?: 0) + 1 }
                    } else {
                        unmatchedIndels.compute(indel) { _, u -> (u ?: 0) + 1 }
                    }
                }
            }
        }

        return result
    }

    private fun SAMRecord.bothEndsInRangeOfCodingTranscripts(): Boolean {
        val thisInRange = codingRegions.any { it.chromosome() == this.contig && this.alignmentStart >= it.start() && this.alignmentStart <= it.end() }
        val mateInRange = codingRegions.any { it.chromosome() == this.contig && this.mateAlignmentStart >= it.start() && this.mateAlignmentStart <= it.end() }

        return thisInRange && mateInRange
    }

}
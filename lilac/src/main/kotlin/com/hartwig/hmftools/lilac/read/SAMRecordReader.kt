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
        const val MAX_DISTANCE = 1000
        val logger = LogManager.getLogger(this::class.java)
        val indelPon = SAMRecordReader::class.java.getResource("/pon/indels.txt")
                .readText()
                .split("\n")
                .map { Indel.fromString(it) }
                .toSet()
    }

    private val codingRegions = transcripts.map { GenomeRegions.create(it.chromosome(), it.codingStart() - MAX_DISTANCE, it.codingEnd() + MAX_DISTANCE) }
    private val unmatchedIndels = mutableMapOf<Indel, Int>()
    private val unmatchedPONIndels = mutableMapOf<Indel, Int>()
    private var alignmentFiltered = 0

    fun alignmentFiltered(): Int {
        return alignmentFiltered
    }

    fun unmatchedIndels(minCount: Int): Map<Indel, Int> {
        return unmatchedIndels.filter { it.value >= minCount }
    }

    fun unmatchedPonIndels(minCount: Int): Map<Indel, Int> {
        return unmatchedIndels.filter { it.value >= minCount }
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
                    val codingRecords = query(variantPosition, codingRegion, bamFile)
                            .filter { recordContainsVariant(variant, it) }

                    val nucleotideFragments = codingRecords
                            .mapNotNull { factory.createAlignmentFragments(it, reverseStrand, codingRegion) }

                    return nucleotideFragments
                }
                hlaCodingRegionOffset += codingRegion.bases().toInt()
            }
        }

        return listOf()
    }

    private fun recordContainsVariant(variant: VariantContextDecorator, record: SAMCodingRecord): Boolean {
        if (variant.alt().length != variant.ref().length) {
            val expectedIndel = Indel(variant.chromosome(), variant.position().toInt(), variant.ref(), variant.alt())
            return record.indels.contains(expectedIndel)
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
        var length = 0
        for (codingRegion in codingRegions) {
            realignedRegions.addAll(realign(codingRegion, reverseStrand, bamFile))
            length += codingRegion.bases().toInt()
        }

        return realignedRegions
                .groupBy { it.id }
                .map { it.value.reduce { x, y -> NucleotideFragment.merge(x, y) } }
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

    private fun query(variantRegion: GenomePosition, nearestCodingRegion: NamedBed, bamFileName: String): List<SAMCodingRecord> {
        val slicer = SAMSlicer(1)
        val result = mutableListOf<SAMCodingRecord>()
        samReaderFactory().open(File(bamFileName)).use { samReader ->
            val consumer = Consumer<SAMRecord> { samRecord ->
                if (samRecord.bothEndsInRangeOfCodingTranscripts()) {
                    result.add(SAMCodingRecord.create(nearestCodingRegion, samRecord))

                } else {
                    alignmentFiltered++
                }
            }

            slicer.slice(variantRegion.chromosome(), variantRegion.position().toInt(), variantRegion.position().toInt(), samReader, consumer)
        }
        return result
    }

    private fun query(codingRegion: NamedBed, bamFileName: String): List<SAMCodingRecord> {
        val slicer = SAMSlicer(1)
        val result = mutableListOf<SAMCodingRecord>()
        samReaderFactory().open(File(bamFileName)).use { samReader ->
            val consumer = Consumer<SAMRecord> { samRecord ->
                if (samRecord.bothEndsInRangeOfCodingTranscripts()) {
                    result.add(SAMCodingRecord.create(codingRegion, samRecord))

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
        for (codingRecord in query(codingRegion, bamFileName)) {
            val fragment = factory.createFragment(codingRecord, reverseStrand, codingRegion)
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
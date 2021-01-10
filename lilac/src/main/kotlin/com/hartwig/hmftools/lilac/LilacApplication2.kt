package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.genome.bed.NamedBed
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.common.genome.region.CodingRegions
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.nuc.AminoAcidCount
import com.hartwig.hmftools.lilac.phase.HeterozygousEvidence
import com.hartwig.hmftools.lilac.phase.PhasedEvidence
import com.hartwig.hmftools.lilac.prot.ProteinSequence
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile.firstFourDigits
import com.hartwig.hmftools.lilac.prot.ProteinSequenceFile.inflate
import com.hartwig.hmftools.lilac.read.Fragment
import com.hartwig.hmftools.lilac.read.SAMRecordRead
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToFirstFourDigits
import org.apache.logging.log4j.LogManager
import java.util.concurrent.Executors
import kotlin.math.min

fun main(args: Array<String>) {
    LilacApplication2().use { x -> x.run() }
}


class LilacApplication2 : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
    }

    private val startTime = System.currentTimeMillis()
    private val transcripts = HmfGenePanelSupplier.allGenesMap37()

    val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("LILAC-%d").build()
    val executorService = Executors.newFixedThreadPool(7, namedThreadFactory)


    val minBaseQual = 30
    val minBaseCount = 2

    override fun run() {


        val minKmerLength = 10
        val maxKmerLength = 20
        val resourcesDir = "/Users/jon/hmf/analysis/hla/resources"
        val sequenceAFile = "${resourcesDir}/A_nuc.txt"
        val sequenceBFile = "${resourcesDir}/B_nuc.txt"
        val sequenceCFile = "${resourcesDir}/C_nuc.txt"
        val bamFile = "/Users/jon/hmf/analysis/hla/GIABvsSELFv004R.hla.bam"
        // 6:29812305-33068473

        val proteinAFile = "${resourcesDir}/A_prot.txt"
        val proteinBFile = "${resourcesDir}/B_prot.txt"
        val proteinCFile = "${resourcesDir}/C_prot.txt"

        LilacApplication.logger.info("Reading protein alignment files")
        val aSequences = ProteinSequenceFile.readFile(proteinAFile).inflate().firstFourDigits()
        val bSequences = ProteinSequenceFile.readFile(proteinBFile).inflate().firstFourDigits()
        val cSequences = ProteinSequenceFile.readFile(proteinCFile).inflate().firstFourDigits()
        val allProteinSequences = mutableListOf<ProteinSequence>()
        allProteinSequences.addAll(aSequences)
        allProteinSequences.addAll(bSequences)
        allProteinSequences.addAll(cSequences)

//        unwrapFile(sequenceAFile)
//        unwrapFile(sequenceBFile)
//        unwrapFile(sequenceCFile)


        // Read from bam
        val reads = mutableListOf<SAMRecordRead>()
        reads.addAll(readFromBam("HLA-A", bamFile))
        reads.addAll(readFromBam("HLA-B", bamFile))
        reads.addAll(readFromBam("HLA-C", bamFile))

        // Join reads together
        val readFragments = reads.groupBy { it.samRecord.readName }.map { Fragment(it.value) }

        // Count individual bases
        val aminoAcidCounts = AminoAcidCount(minBaseQual, readFragments)
        val excludedIndices = setOf(24, 114, 206, 298, 337, 348, 349, 362, 364, 365, 366)
        val heterozygousIndices = aminoAcidCounts.heterozygousIndices(minBaseCount).filter { it !in excludedIndices }
        val initialCandidates = initialCandidates(excludedIndices, aminoAcidCounts, allProteinSequences)
        println("${allProteinSequences.size} types")
        println("${initialCandidates.size} candidates after amino acid filtering")

        val heterozygousEvidence = HeterozygousEvidence(minBaseQual,
                heterozygousIndices.filter { it !in (298..370) }, readFragments)

        val initialEvidence = heterozygousEvidence.initialEvidence()


        println("HERE")
        println(PhasedEvidence.evidence(minBaseQual, readFragments, 93, 94, 96))
//        println(PhasedEvidence2.evidence(minBaseQual, readFragments, 90, 92, 93, 94, 96, 99, 100, 102, 103))
//        println(PhasedEvidence2.evidence(minBaseQual, readFragments, 67, 68, 75, 85, 86, 88, 89, 90))
//        println(PhasedEvidence2.evidence(minBaseQual, readFragments, 67, 68, 75, 85, 86, 88, 89, 90, 92, 93, 94, 96, 99, 100, 102, 103))
//

        var candidates = initialCandidates

        val allEvidence = mutableSetOf<PhasedEvidence>()
        allEvidence.addAll(initialEvidence)

        var unprocessedEvidence = initialEvidence

        for (i in 0..1) {
            if (unprocessedEvidence.isEmpty()) {
                break
            }

            val topEvidence = unprocessedEvidence[0]

            candidates = matchingCandidates(topEvidence, candidates)

            if (i % 1 == 0) {
                println("Iteration ${i}: ${candidates.size} candidates includes ${checkCandidates(candidates)} actual using evidence with ${topEvidence.evidence.size} types -> $topEvidence ")
            }

            val newEvidence = heterozygousEvidence.extendEvidence(topEvidence, allEvidence)
            allEvidence.addAll(newEvidence)

            val updatedEvidence = mutableSetOf<PhasedEvidence>()
            updatedEvidence.addAll(unprocessedEvidence.drop(1))
            updatedEvidence.addAll(newEvidence)

            unprocessedEvidence = updatedEvidence.sorted()


        }



        println(heterozygousIndices)



        println("Fragments overlapping")
//        for (heterozygousIndex in heterozygousIndices) {
//            var intersectionCount = 0
//            for (fragment in readFragments) {
//                if (fragment.aminoAcidIndices().contains(heterozygousIndex)) {
//                    intersectionCount += (fragment.aminoAcidIndices() intersect heterozygousIndices).size
//                }
//            }
//
//            println("Insersection count $heterozygousIndex -> $intersectionCount")
//
//        }

//        for (fragment in readFragments) {
//            if (fragment.aminoAcidIndices().contains(1)) {
//                println(fragment.aminoAcidIndices() intersect heterozygousIndices)
//            }
//        }


//        println("Phased Evidence")
//        val jon2 = PhasedEvidence2.evidence(minBaseQual, readFragments, 85, 86, 88, 89, 90, 92)
////        val jon2 = PhasedEvidence2.evidence(minBaseQual, readFragments, 200, 203, 205, 207, 212, 217)
//        val jon2Candidates = matchingCandidates(jon2, initialCandidates)
//        println("${jon2Candidates.size} candidates after p[ahsed filtering")
//        checkCandidates(jon2Candidates)
//
//        println(jon2)


//        val jonjon = reads.groupBy { it.samRecord.readName }.filter { it.value.size > 1 }
//        println(jonjon.size)
//
//
//        var jon = PhasedSequence.locallyPhased(0, 23, reads)
//        for (phasedSequence in jon) {
//            println(phasedSequence)
//        }
//
//        jon = PhasedSequence.locallyPhased(7, 10, reads)
//        for (phasedSequence in jon) {
//            println(phasedSequence)
//        }
//
//
//        var aaCandidates = allProteinSequences.toSet()
//
//        for (i in 0..360) {
//            aaCandidates = reduceAACandidates(i, aaCandidates, aminoAcidCounts)
//        }
//        val exon1 = setOf("MAVMAPRTLLLLLSGALALTQTWA", "MRVMAPRTLILLLSGALALTETWA", "MRVMAPRALLLLLSGGLALTETWA", "MRVTAPRTLLLLLWGALALTETWA", "MLVMAPRTVLLLLSAALALTETWA")
//
//        println(aaCandidates.size)
//        aaCandidates = aaCandidates.filter { it.fullProteinSequence.substring(0, 24) in exon1 }.toSet()
//        println(aaCandidates.size)
//
////
////
//        println("A: -> " + aaCandidates.filter { it.allele.toString().startsWith("A") }.count())
//        println("B: -> " + aaCandidates.filter { it.allele.toString().startsWith("B") }.count())
//        println("C: -> " + aaCandidates.filter { it.allele.toString().startsWith("C") }.count())
////
//
//        aminoAcidCount.writeVertically("/Users/jon/hmf/analysis/hla/aminoAcidCount.qual30.txt")
//


//        val specificSequencesA = setOf(HlaAllele("A*01:01:01:01"), HlaAllele("A*11:01:01:01"))
//        val specificSequencesB = setOf(HlaAllele("B*08:01:01:01"), HlaAllele("B*56:01:01:01"))
//        val specificSequencesC = setOf(HlaAllele("C*01:02:01:01"), HlaAllele("C*07:01:01:01"))
//
//        val nucleotideCount = NucleotideCount(1101)
//        for (realignedRegion in realignedRegions) {
//            for (i in realignedRegion.nucleotideIndices) {
//
//                when (realignedRegion.charAt(i, minBaseQual)) {
//                    'G' -> nucleotideCount.gCount[i] += 1
//                    'A' -> nucleotideCount.aCount[i] += 1
//                    'T' -> nucleotideCount.tCount[i] += 1
//                    'C' -> nucleotideCount.cCount[i] += 1
//                }
//            }
//        }
//        nucleotideCount.writeVertically("/Users/jon/hmf/analysis/hla/nucleotideCount.qual30.txt")
//
//        println("HERE")

    }


    private fun checkCandidates(candidates: Collection<ProteinSequence>): Int {
        var count = 0

        if (candidates.any { it.allele == HlaAllele("A*01:01:01:01") }) {
            count++;
        }

        if (candidates.any { it.allele == HlaAllele("A*11:01:01:01") }) {
            count++;
        }
        if (candidates.any { it.allele == HlaAllele("B*08:01:01:01") }) {
            count++;
        }
        if (candidates.any { it.allele == HlaAllele("B*56:01:01:01") }) {
            count++;
        }
        if (candidates.any { it.allele == HlaAllele("C*01:02:01:01") }) {
            count++;
        }
        if (candidates.any { it.allele == HlaAllele("C*07:01:01:01") }) {
            count++;
        }


        return count;
    }

    private fun initialCandidates(excludedLocations: Collection<Int>, aminoAcidCount: AminoAcidCount, candidates: List<ProteinSequence>): List<ProteinSequence> {
        var result = candidates
        val locations = (0 until min(360, aminoAcidCount.length)).toSet() subtract excludedLocations
        for (location in locations) {
            result = filterCandidates(location, aminoAcidCount.aminoAcidAt(location, minBaseCount), result)
        }
        return result
    }

    private fun filterCandidates(index: Int, aminoAcids: Collection<Char>, candidates: Collection<ProteinSequence>): List<ProteinSequence> {
        return candidates.filter { it.length <= index || it.fullProteinSequence[index] == '*' || it.fullProteinSequence[index] in aminoAcids }
    }

    private fun matchingCandidates(evidence: PhasedEvidence, candidates: Collection<ProteinSequence>): List<ProteinSequence> {
        return candidates.filter { it.consistentWith(evidence) }
    }


    private fun readFromBam(gene: String, bamFile: String): List<SAMRecordRead> {
        val transcript = transcripts[gene]!!
        val reverseStrand = transcript.strand() == Strand.REVERSE
        val codingRegions = if (reverseStrand) codingRegions(gene).reversed() else codingRegions(gene)


        println("Processing ${transcript.gene()} (${transcript.chromosome()}:${transcript.codingStart()}-${transcript.codingEnd()})")
        val realignedRegions = mutableListOf<SAMRecordRead>()
        var length = 0
        for (codingRegion in codingRegions) {
            realignedRegions.addAll(SAMRecordRead.realign(length, codingRegion, reverseStrand, bamFile))
            length += codingRegion.bases().toInt()
//            println((length - 1) / 3)
        }
        return realignedRegions
    }

    fun codingRegions(gene: String): List<NamedBed> {
        return CodingRegions.codingRegions(transcripts[gene]!!)
    }


    fun unwrapFile(fileName: String) {
        val input = HlaSequenceFile.readFile(fileName)
        val output = input.reduceToFirstFourDigits()
        val inflated = output.map { it.inflate(input[0].rawSequence) }
        val deflated = inflated.map { it.deflate(inflated[0].rawSequence) }
        HlaSequenceFile.writeFile(fileName.replace(".txt", ".unwrapped.txt"), output)
        HlaSequenceFile.writeFile(fileName.replace(".txt", ".deflated.txt"), deflated)
    }


    override fun close() {
        executorService.shutdown()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }

}
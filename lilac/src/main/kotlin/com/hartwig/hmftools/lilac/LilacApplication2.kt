package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.hla.HlaAlleleCount
import com.hartwig.hmftools.lilac.nuc.SequenceCount
import com.hartwig.hmftools.lilac.phase.HeterozygousEvidence
import com.hartwig.hmftools.lilac.phase.PhasedEvidence
import com.hartwig.hmftools.lilac.read.Fragment
import com.hartwig.hmftools.lilac.read.FragmentSequences
import com.hartwig.hmftools.lilac.read.FragmentSequencesFile
import com.hartwig.hmftools.lilac.read.SAMRecordRead
import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.deflate
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.inflate
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


    private fun filterCandidatesOnNucleotides(minEvidence: Int, candidates: Collection<HlaSequence>, fragments: List<Fragment>, vararg nucleotides: Int): List<HlaSequence> {
        try {
            val reads = fragments
                    .filter { it.containsAllNucleotides(*nucleotides) }
                    .map { it.nucleotides(*nucleotides) }
                    .groupingBy { it }
                    .eachCount()
                    .filter { it.value >= minEvidence }
                    .keys
                    .map { it.toCharArray() }


            return candidates.filter { it.consistentWith(nucleotides, reads) }
        } catch (e: Exception) {
            throw e
        }
    }

    private fun filterCandidatesOnExonBoundaryNucleotide(aminoAcidIndex: Int, minEvidence: Int, candidates: Collection<HlaSequence>, fragments: List<Fragment>): List<HlaSequence> {
        val firstBaseCandidates = filterCandidatesOnNucleotides(minEvidence, candidates, fragments, aminoAcidIndex * 3)
        return filterCandidatesOnNucleotides(minEvidence, firstBaseCandidates, fragments, aminoAcidIndex * 3 + 1, aminoAcidIndex * 3 + 2)
    }

    private fun filterCandidatesOnExonBoundaries(exonBoundaries: Collection<Int>, minEvidence: Int, candidates: Collection<HlaSequence>, fragments: List<Fragment>): List<HlaSequence> {
        var result = candidates.toList()
        for (exonBoundary in exonBoundaries) {
            result = filterCandidatesOnExonBoundaryNucleotide(exonBoundary, minEvidence, result, fragments)
        }

        return result
    }


    override fun run() {
        logger.info("Starting")

        val minKmerLength = 10
        val maxKmerLength = 20
        val resourcesDir = "/Users/jon/hmf/analysis/hla/resources"
        val bamFile = "/Users/jon/hmf/analysis/hla/GIABvsSELFv004R.hla.bam"


        logger.info("Reading bam $bamFile")
        val readFragments = readFromBam(bamFile)


        // Count individual bases
        val aminoAcidCounts = SequenceCount.aminoAcids(minBaseQual, readFragments)
        aminoAcidCounts.writeVertically("/Users/jon/hmf/analysis/hla/aminoacids.count.txt")

        val nucleotideCounts = SequenceCount.nucleotides(minBaseQual, readFragments)
        nucleotideCounts.writeVertically("/Users/jon/hmf/analysis/hla/nucleotides.count.txt")

        val aminoAcidBoundaries = setOf(24, 114, 206, 298, 337, 348, 349, 362, 364, 365, 366)

        val aBoundaries = setOf(24, 114, 206, 298, 337, 348, 364, 365)
        val bBoundaries = setOf(24, 114, 206, 298, 337, 348, 362)
        val cBoundaries = setOf(24, 114, 206, 298, 338, 349, 365, 366)
        val allBoundaries = (aBoundaries + bBoundaries + cBoundaries)//.filter { it !in setOf(24, 114) }
        val commonBoundaries = aBoundaries intersect bBoundaries intersect cBoundaries

        val excludedIndices = allBoundaries.toSet()

        val hetLoci = aminoAcidCounts.heterozygousIndices(minBaseCount)
        val heterozygousIndices = hetLoci
                .filter { it !in excludedIndices }

        println("Heterozygous locations")
        println(heterozygousIndices)

        LilacApplication.logger.info("Reading nucleotide files")
        val nucleotideSequences = readSequenceFiles { "${resourcesDir}/${it}_nuc.txt" }

        logger.info("Reading protein files")
        val allProteinSequences = readSequenceFiles { "${resourcesDir}/${it}_prot.txt" }

//        val initialNucleotideCandidates = initialNucleotideCandidates(nucleotideCounts, nucleotideSequences)
        val boundaryNucleotideCandidates = filterCandidatesOnExonBoundaries(allBoundaries, 1, nucleotideSequences, readFragments)
        val nucleotideFilteredAlleles = boundaryNucleotideCandidates.map { it.contig }
        println("${nucleotideSequences.size} types")
        println("${boundaryNucleotideCandidates.size} candidates after nucleotide filtering")


        val initialCandidates = initialCandidates(excludedIndices, aminoAcidCounts, allProteinSequences)
                .filter { it.contig in nucleotideFilteredAlleles }
        println("${allProteinSequences.size} types")
        println("${initialCandidates.size} candidates after amino acid filtering")

        // Specific Evidence
        println(PhasedEvidence.evidence(30, readFragments, 32, 67))

        // Full Evidence
        val consecutiveEvidence = consecutiveEvidence(excludedIndices, aminoAcidCounts, readFragments, initialCandidates)
        println("Constructed ${consecutiveEvidence.size} consecutive sequences")

        var candidates = initialCandidates
        for (i in consecutiveEvidence.indices) {
            val evidence = consecutiveEvidence[i]
            candidates = matchingCandidates(evidence, candidates)
            println("$i ->  ${candidates.size} candidates includes ${checkCandidates(candidates)} actual types -> $evidence ")
        }

        val fragmentSequences = FragmentSequences.create(readFragments, hetLoci, candidates)

        FragmentSequencesFile.writeFile("/Users/jon/hmf/analysis/hla/fragments.txt", fragmentSequences)

        val groupCoverage = HlaAlleleCount.groupCoverage(fragmentSequences)
        val proteinCoverage = HlaAlleleCount.proteinCoverage(fragmentSequences)
        println("SDFSD")

        println("Group Coverage")
        for (hlaAlleleCount in groupCoverage) {
            println(hlaAlleleCount)
        }

        println("Protein Coverage")
        for (hlaAlleleCount in proteinCoverage) {
            println(hlaAlleleCount)
        }

//        val jon = PhasedEvidence.combineOverlapping(consecutiveEvidence[2], consecutiveEvidence[3])
//        val exon2 = PhasedEvidence.combineOverlapping(jon, consecutiveEvidence[4])
//
//        println(consecutiveEvidence[2])
//        println(consecutiveEvidence[3])
//        println(consecutiveEvidence[4])
//        println(jon)
//        println(exon2)
//
//        println("FILTERINUG")
//
//
//        var sparseEvidence = PhasedEvidence.evidence(30, readFragments, 29, 32, 67)
//        println(sparseEvidence)
//        var jon3 = exon2.reduceInLightOfNewEvidence(sparseEvidence)
//        println(jon3)
////
////
//        sparseEvidence = PhasedEvidence.evidence(30, readFragments, 32, 34, 67)
//        println(sparseEvidence)
//        jon3 = jon3.reduceInLightOfNewEvidence(sparseEvidence)
//        println(jon3)
//
//        sparseEvidence = PhasedEvidence.evidence(30, readFragments, 29,89)
//        println(sparseEvidence)
//        jon3 = jon3.reduceInLightOfNewEvidence(sparseEvidence)
//        println(jon3)
//
//        sparseEvidence = PhasedEvidence.evidence(30, readFragments, 32,35,68)
//        println(sparseEvidence)
//        jon3 = jon3.reduceInLightOfNewEvidence(sparseEvidence)
//        println(jon3)
//
//        sparseEvidence = PhasedEvidence.evidence(30, readFragments, 47,68, 90)
//        println(sparseEvidence)
//        jon3 = jon3.reduceInLightOfNewEvidence(sparseEvidence)
//        println(jon3)


//        var exon2Pairs = mutableListOf<PhasedEvidence>()
//        var jon3 = exon2
//        for (i in exon2.aminoAcidIndices) {
//            for (j in exon2.aminoAcidIndices) {
////                for (k in exon2.aminoAcidIndices) {
//
//                    if (j > i) {
//                        var newEvidence = PhasedEvidence.evidence(30, readFragments, i, j)
//                        if (newEvidence.totalEvidence() > 20) {
//                            println(newEvidence)
//                            jon3 = jon3.reduceInLightOfNewEvidence(newEvidence)
//                            println(jon3)
//                        }
//                        exon2Pairs.add(newEvidence)
////                    }
//                }
//            }
//        }
//        exon2Pairs.sort()


//        println("${candidates.size} candidates after full match filtering")


//        println("AFTER PARTIAL MATCHING")
//
//        var consolidated = consolidate(fullEvidence)
//        consolidated = consolidate(consolidated)
//        consolidated = consolidate(consolidated)
//        consolidated = consolidate(consolidated)
//        consolidated = consolidate(consolidated)
//        consolidated = longestFullEvidence(consolidated)
//
//        for (i in consolidated.indices) {
//            val evidence = consolidated[i]
//            candidates = matchingCandidates(evidence, candidates)
//            println("$i ->  ${candidates.size} candidates includes ${checkCandidates(candidates)} actual types -> $evidence ")
//        }
//        println("${candidates.size} candidates after partial match filtering")


        val sequences = candidates.map { HlaSequence(it.contig, it.sequence) }
        HlaSequenceFile.writeFile("/Users/jon/hmf/analysis/hla/candidates.inflate.txt", sequences)
        HlaSequenceFile.wipeFile("/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.writeBoundary(aBoundaries, "/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.writeBoundary(bBoundaries, "/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.writeBoundary(cBoundaries, "/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.appendFile("/Users/jon/hmf/analysis/hla/candidates.deflate.txt", sequences.deflate())


//        var combined = PhasedEvidence.combineOverlapping(fullEvidence[19], fullEvidence[20])
//        println(combined)
//        combined = PhasedEvidence.combineOverlapping(combined, fullEvidence[21])
//        println(combined)
//        combined = PhasedEvidence.combineOverlapping(combined, fullEvidence[22])
//        println(combined)
//        println("SDF")
//
//        combined = PhasedEvidence.combineOverlapping(fullEvidence[22], fullEvidence[23])
//        println(combined)


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

    private fun consolidate(evidence: List<PhasedEvidence>): List<PhasedEvidence> {
        val result = mutableListOf<PhasedEvidence>()
        val used = mutableSetOf<Int>()

        for (i in evidence.indices) {
            for (j in i + 1 until evidence.size) {
                if (evidence[i].overlaps(evidence[j])) {
                    val consolidated = PhasedEvidence.combineOverlapping(evidence[i], evidence[j])
                    if (consolidated.evidence.size <= 6) {
                        result.add(consolidated)
                        used.add(i)
                        used.add(j)
                    }
                }
            }

        }

//        for (i in evidence.indices subtract used) {
//            result.add(evidence[i])
//        }

        result.addAll(evidence)

        return result.toSet().toList().sortedBy { it.aminoAcidIndices[0] }
    }
//

    private fun consecutiveEvidence(excludedIndices: Set<Int>, aminoAcidCounts: SequenceCount, readFragments: List<Fragment>, initialCandidates: List<HlaSequence>): List<PhasedEvidence> {

        var candidates = initialCandidates
        val heterozygousIndices = aminoAcidCounts.heterozygousIndices(minBaseCount).filter { it !in excludedIndices }
        val heterozygousEvidence = HeterozygousEvidence(minBaseQual, heterozygousIndices, readFragments)

        val allEvidence = mutableSetOf<PhasedEvidence>()
        val initialEvidence = heterozygousEvidence.consecutiveEvidence()
        var unprocessedEvidence = initialEvidence

        allEvidence.addAll(initialEvidence)

        var i = 0
        while (unprocessedEvidence.isNotEmpty()) {
            val topEvidence = unprocessedEvidence[0]
            allEvidence.add(topEvidence)

            candidates = matchingCandidates(topEvidence, candidates)
//            println("${i++} ->  ${candidates.size} candidates includes ${checkCandidates(candidates)} actual types -> $topEvidence")


            val newEvidence = heterozygousEvidence.extendConsecutive(topEvidence, allEvidence)
            allEvidence.addAll(newEvidence)

            val updatedEvidence = mutableSetOf<PhasedEvidence>()
            updatedEvidence.addAll(unprocessedEvidence
                    .drop(1)
//                    .filter { !newEvidence.any { x -> x.contains(it) } }
            )
            updatedEvidence.addAll(newEvidence)

            unprocessedEvidence = updatedEvidence.sorted()
        }

        return longestFullEvidence(allEvidence)

    }

    private fun longestFullEvidence(evidence: Collection<PhasedEvidence>): List<PhasedEvidence> {

        fun Collection<PhasedEvidence>.otherContains(victim: PhasedEvidence): Boolean {
            return this.any { it != victim && it.contains(victim) }
        }

        return evidence
                .filter { !evidence.otherContains(it) }
                .sortedBy { it.aminoAcidIndices[0] }
    }


    private fun checkCandidates(candidates: Collection<HlaSequence>): Int {
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

    private fun initialCandidates(excludedLocations: Collection<Int>, aminoAcidCount: SequenceCount, candidates: List<HlaSequence>): List<HlaSequence> {
        var result = candidates
        val locations = (0 until min(360, aminoAcidCount.length)).toSet() subtract excludedLocations
        for (location in locations) {
            result = filterCandidates(location, aminoAcidCount.sequenceAt(location, minBaseCount), result)
        }
        return result
    }

    private fun initialNucleotideCandidates(aminoAcidCount: SequenceCount, candidates: List<HlaSequence>): List<HlaSequence> {
        var result = candidates
        val locations = (0 until min(1080, aminoAcidCount.length)).toSet()
        for (location in locations) {
            result = filterCandidates(location, aminoAcidCount.sequenceAt(location, minBaseCount), result)
        }
        return result
    }

    private fun filterCandidates(index: Int, expectedCharacters: Collection<Char>, candidates: Collection<HlaSequence>): List<HlaSequence> {
        return candidates.filter { it.length <= index || it.sequence[index] == '*' || it.sequence[index] in expectedCharacters }
    }

    private fun matchingCandidates(evidence: PhasedEvidence, candidates: Collection<HlaSequence>): List<HlaSequence> {
        return candidates.filter { it.consistentWith(evidence) }
    }

    private fun readFromBam(bamFile: String): List<Fragment> {
        val reads = mutableListOf<SAMRecordRead>()
        reads.addAll(SAMRecordRead.readFromBam(transcripts["HLA-A"]!!, bamFile))
        reads.addAll(SAMRecordRead.readFromBam(transcripts["HLA-B"]!!, bamFile))
        reads.addAll(SAMRecordRead.readFromBam(transcripts["HLA-C"]!!, bamFile))

        return Fragment.fromReads(minBaseQual, reads)
    }


    fun unwrapFile(fileName: String) {
        val input = HlaSequenceFile.readFile(fileName)
        val output = input.reduceToFirstFourDigits()
        val inflated = output.map { it.inflate(input[0].rawSequence) }
        val deflated = inflated.map { it.deflate(inflated[0].rawSequence) }
        HlaSequenceFile.writeFile(fileName.replace(".txt", ".unwrapped.txt"), output)
        HlaSequenceFile.writeFile(fileName.replace(".txt", ".deflated.txt"), deflated)
    }

    private fun readSequenceFiles(filenameSupplier: (Char) -> String): List<HlaSequence> {

        val aFile = filenameSupplier('A')
        val bFile = filenameSupplier('B')
        val cFile = filenameSupplier('C')

        val aSequence = HlaSequenceFile.readFile(aFile).inflate().reduceToFirstFourDigits()
        val bSequence = HlaSequenceFile.readFile(bFile).inflate().reduceToFirstFourDigits()
        val cSequence = HlaSequenceFile.readFile(cFile).inflate().reduceToFirstFourDigits()

        val result = mutableListOf<HlaSequence>()
        result.addAll(aSequence)
        result.addAll(bSequence)
        result.addAll(cSequence)

        val maxLength = result.map { it.sequence.length }.max()!!

        return result
                .filter { it.sequence.isNotEmpty() }
//                .map { it.pad(maxLength) }
    }


    override fun close() {
        executorService.shutdown()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }

}
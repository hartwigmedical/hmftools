package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.hla.HlaAlleleCoverageFactory
import com.hartwig.hmftools.lilac.hla.HlaAlleleCoverageFactory.Companion.coverageString
import com.hartwig.hmftools.lilac.hla.HlaComplex
import com.hartwig.hmftools.lilac.nuc.SequenceCount
import com.hartwig.hmftools.lilac.phase.ExtendedEvidence
import com.hartwig.hmftools.lilac.phase.PhasedEvidence
import com.hartwig.hmftools.lilac.phase.TypeEvidence
import com.hartwig.hmftools.lilac.read.*
import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.deflate
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.inflate
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.specificProteins
import org.apache.logging.log4j.LogManager
import java.util.concurrent.Executors

fun main(args: Array<String>) {
    LilacApplication2().use { x -> x.run() }
}


class LilacApplication2 : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
        const val HLA_A = "HLA-A"
        const val HLA_B = "HLA-B"
        const val HLA_C = "HLA-C"
    }

    private val startTime = System.currentTimeMillis()
    private val transcripts = HmfGenePanelSupplier.allGenesMap37()

    val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("LILAC-%d").build()
    val executorService = Executors.newFixedThreadPool(7, namedThreadFactory)


    val minBaseQual = 30
    val minBaseCount = 2


    override fun run() {
        logger.info("Starting")

        val aProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348, 364, 365)
        val bProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348, 362)
        val cProteinExonBoundaries = setOf(24, 114, 206, 298, 338, 349, 365, 366)
        val allProteinExonBoundaries = (aProteinExonBoundaries + bProteinExonBoundaries + cProteinExonBoundaries)
        val allNucleotideExonBoundaryStarts = allProteinExonBoundaries.map { it * 3 }



        val resourcesDir = "/Users/jon/hmf/analysis/hla/resources"
        val bamFile = "/Users/jon/hmf/analysis/hla/GIABvsSELFv004R.hla.bam"


        logger.info("Reading nucleotides from  $bamFile")
        val rawNucleotideFragments = readFromBam(bamFile)
        val nucleotideFragmentFactory = NucleotideFragmentFactory(minBaseCount, rawNucleotideFragments, aProteinExonBoundaries, bProteinExonBoundaries, cProteinExonBoundaries)

        val nucleotideCounts = SequenceCount.nucleotides(minBaseCount, rawNucleotideFragments)
        val nucleotideHeterozygousLoci = nucleotideCounts.heterozygousIndices()
        val nucleotideFragments = NucleotideFragmentEnrichment(allNucleotideExonBoundaryStarts, nucleotideCounts)
                .enrichHomSpliceJunctions(rawNucleotideFragments)


        val aminoAcidFragments = nucleotideFragments.map { it.toAminoAcidFragment() }
        val aminoAcidCounts = SequenceCount.aminoAcids(minBaseCount, aminoAcidFragments)




        aminoAcidCounts.writeVertically("/Users/jon/hmf/analysis/hla/aminoacids.count.txt")
        nucleotideCounts.writeVertically("/Users/jon/hmf/analysis/hla/nucleotides.count.txt")


//.filter { it !in setOf(24, 114) }
        val nucleotideLoci = allProteinExonBoundaries.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) } intersect nucleotideHeterozygousLoci
        logger.info("Het nucleotide exon boundaries: " + nucleotideLoci)

        val commonBoundaries = aProteinExonBoundaries intersect bProteinExonBoundaries intersect cProteinExonBoundaries

        val excludedIndices = allProteinExonBoundaries.toSet()

        val hetLoci = aminoAcidCounts.heterozygousIndices()
        val heterozygousIndices = hetLoci
                .filter { it !in excludedIndices }

        println("Heterozygous locations")
        println(heterozygousIndices)

        LilacApplication.logger.info("Reading nucleotide files")
        val nucleotideSequences = readNucleotideFiles(resourcesDir)

        logger.info("Reading protein files")
        val allProteinSequences = readProteinFiles(resourcesDir)

        val boundaryNucleotideCandidates = filterCandidatesOnExonBoundaries(allProteinExonBoundaries, 2, nucleotideSequences, aminoAcidFragments)
        val nucleotideFilteredAlleles = boundaryNucleotideCandidates.map { it.contig }
        println("${nucleotideSequences.size} types")
        println("${boundaryNucleotideCandidates.size} candidates after nucleotide filtering")


        val initialCandidates = initialCandidates(excludedIndices, aminoAcidCounts, allProteinSequences)
                .filter { it.contig in nucleotideFilteredAlleles }
        println("${allProteinSequences.size} types")
        println("${initialCandidates.size} candidates after amino acid filtering")

        val typeEvidenceFactory = TypeEvidence(minBaseCount, 20, rawNucleotideFragments, aProteinExonBoundaries, bProteinExonBoundaries, cProteinExonBoundaries)

        val typeAEvidence = typeEvidenceFactory.typeAEvidence()
        val typeBEvidence = typeEvidenceFactory.typeBEvidence()
        val typeCEvidence = typeEvidenceFactory.typeCEvidence()


        val aCandidates = filterCandidates(initialCandidates.filter { it.allele.gene == "A" }, typeAEvidence)
        val bCandidates = filterCandidates(initialCandidates.filter { it.allele.gene == "B" }, typeBEvidence)
        val cCandidates = filterCandidates(initialCandidates.filter { it.allele.gene == "C" }, typeCEvidence)
        val consecutiveEvidenceCandidates = aCandidates + bCandidates + cCandidates

        for (consecutiveEvidenceCandidate in consecutiveEvidenceCandidates) {
            println(consecutiveEvidenceCandidate)
        }

        val candidateAlleles = consecutiveEvidenceCandidates.map { it.allele }
        val candidateAlleleSpecificProteins = candidateAlleles.map { it.specificProtein() }

        val aminoAcidCandidates = allProteinSequences.filter { it.allele in candidateAlleles }
        val nucleotideCandidates = nucleotideSequences.filter { it.allele.specificProtein() in candidateAlleleSpecificProteins }

        val coverageFactory = HlaAlleleCoverageFactory(aminoAcidFragments, hetLoci, aminoAcidCandidates, nucleotideLoci, nucleotideCandidates)
        val proteinCoverage = coverageFactory.proteinCoverage(candidateAlleles)


//        val fragmentSequences = FragmentAlleles.create(readFragments, hetLoci, consecutiveEvidenceCandidates)


//        FragmentSequencesFile.writeFile("/Users/jon/hmf/analysis/hla/fragments.txt", fragmentSequences)

        val groupCoverage = coverageFactory.groupCoverage(candidateAlleles)
        println("SDFSD")

        val confirmedGroups = listOf(HlaAllele("B*08"), HlaAllele("C*07"), HlaAllele("C*01"), HlaAllele("B*56"), HlaAllele("A*11"), HlaAllele("A*01"))
        val confimedProtein = listOf(HlaAllele("C*01:02:01:01"), HlaAllele("A*01:01:01:01"), HlaAllele("B*56:01:01:01"))
        val complexes = HlaComplex.complexes(confirmedGroups, confimedProtein, consecutiveEvidenceCandidates.map { it.allele })
        println("TotalCoverage\tUniqueCoverage\tSharedCoverage\tAllele1\tAllele2\tAllele3\tAllele4\tAllele5\tAllele6")
        for (complex in complexes) {
            val complexCoverage = coverageFactory.proteinCoverage(complex.alleles)
            println(complexCoverage.coverageString(confimedProtein))
        }


        val testCoverage = coverageFactory.proteinCoverage(
                listOf(HlaAllele("C*01:02:01:01"), HlaAllele("A*01:01:01:01"), HlaAllele("B*56:01:01:01"), HlaAllele("B*08:01:01:01"), HlaAllele("C*07:01:01:01"),

                        HlaAllele("A*11:01:01:01"), HlaAllele("A*11:303"), HlaAllele("A*11:353")))
        println(testCoverage)


        val sequences = consecutiveEvidenceCandidates.map { HlaSequence(it.contig, it.sequence) }
        HlaSequenceFile.writeFile("/Users/jon/hmf/analysis/hla/candidates.inflate.txt", sequences)
        HlaSequenceFile.wipeFile("/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.writeBoundary(aProteinExonBoundaries, "/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.writeBoundary(bProteinExonBoundaries, "/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.writeBoundary(cProteinExonBoundaries, "/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.appendFile("/Users/jon/hmf/analysis/hla/candidates.deflate.txt", sequences.deflate())


    }
//

    private fun filterCandidates(initialCandidates: List<HlaSequence>, evidence: List<PhasedEvidence>): List<HlaSequence> {
        var candidates = initialCandidates
        for (i in evidence.indices) {
            val newEvidence = evidence[i]
            candidates = matchingCandidates(newEvidence, candidates)
//            println("$i ->  ${candidates.size} candidates includes ${checkCandidates(candidates)} actual types -> $newEvidence ")
        }

        return candidates
    }


    private fun consecutiveEvidence(
            aminoAcidBoundaries: Set<Int>,
            aminoAcidCountsOld: SequenceCount,
            rawNucleotideFragments: List<NucleotideFragment>,
            aExcludedAminoAcidReads: Set<Int>,
            bExcludedAminoAcidReads: Set<Int>,
            cExcludedAminoAcidReads: Set<Int>,
            initialCandidates: List<HlaSequence>): List<PhasedEvidence> {

        fun exclude(fragment: NucleotideFragment, gene: String, excludedNucleotides: Collection<Int>): Boolean {
            return fragment.alignedGene == gene && excludedNucleotides.any { fragment.containsNucleotide(it) }
        }

        val aExcludedNucleotides = aExcludedAminoAcidReads.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }
        val bExcludedNucleotides = bExcludedAminoAcidReads.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }
        val cExcludedNucleotides = cExcludedAminoAcidReads.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }


        val filteredNucleotideFragments = rawNucleotideFragments
                .filter { !exclude(it, HLA_A, aExcludedNucleotides) && !exclude(it, HLA_B, bExcludedNucleotides) && !exclude(it, HLA_C, cExcludedNucleotides) }
        val filteredNucleotideFragmentCounts = SequenceCount.nucleotides(minBaseCount, filteredNucleotideFragments)


        println("ENRICHING")
        val enrichedNucleotideFragments = NucleotideFragmentEnrichment(aminoAcidBoundaries.map { it * 3 }, filteredNucleotideFragmentCounts).enrichHomSpliceJunctions(filteredNucleotideFragments)
        val aminoAcidFragments = enrichedNucleotideFragments.map { it.toAminoAcidFragment() }
        val aminoAcidCounts = SequenceCount.aminoAcids(minBaseCount, aminoAcidFragments)

         var candidates = initialCandidates
        val heterozygousIndices = aminoAcidCounts.heterozygousIndices()//.filter { it !in aminoAcidBoundaries }
        val heterozygousEvidence = ExtendedEvidence(2, 20, heterozygousIndices, aminoAcidFragments)

        val allEvidence = mutableSetOf<PhasedEvidence>()
        val initialEvidence = heterozygousEvidence.initialEvidence()
        var unprocessedEvidence = initialEvidence

        allEvidence.addAll(initialEvidence)

        val jon = PhasedEvidence.evidence(aminoAcidFragments, 337, 338)

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
        val locations = (0 until aminoAcidCount.length).toSet() subtract excludedLocations
        for (location in locations) {
            result = filterCandidates(location, aminoAcidCount.sequenceAt(location), result)
        }
        return result
    }

    private fun initialNucleotideCandidates(aminoAcidCount: SequenceCount, candidates: List<HlaSequence>): List<HlaSequence> {
        var result = candidates
        val locations = (0 until aminoAcidCount.length).toSet()
        for (location in locations) {
            result = filterCandidates(location, aminoAcidCount.sequenceAt(location), result)
        }
        return result
    }

    private fun filterCandidates(index: Int, expectedCharacters: Collection<Char>, candidates: Collection<HlaSequence>): List<HlaSequence> {
        return candidates.filter { it.length <= index || it.sequence[index] == '*' || it.sequence[index] in expectedCharacters }
    }

    private fun matchingCandidates(evidence: PhasedEvidence, candidates: Collection<HlaSequence>): List<HlaSequence> {
        return candidates.filter { it.consistentWith(evidence) }
    }

    private fun readFromBam(bamFile: String): List<NucleotideFragment> {
        val reads = mutableListOf<SAMRecordRead>()
        reads.addAll(SAMRecordRead.readFromBam(transcripts[HLA_A]!!, bamFile))
        reads.addAll(SAMRecordRead.readFromBam(transcripts[HLA_B]!!, bamFile))
        reads.addAll(SAMRecordRead.readFromBam(transcripts[HLA_C]!!, bamFile))

        return NucleotideFragment.fromReads(minBaseQual, reads)
    }


    fun unwrapFile(fileName: String) {
        val input = HlaSequenceFile.readFile(fileName)
        val output = input.specificProteins()
        val inflated = output.map { it.inflate(input[0].rawSequence) }
        val deflated = inflated.map { it.deflate(inflated[0].rawSequence) }
        HlaSequenceFile.writeFile(fileName.replace(".txt", ".unwrapped.txt"), output)
        HlaSequenceFile.writeFile(fileName.replace(".txt", ".deflated.txt"), deflated)
    }

    private fun readNucleotideFiles(resourcesDir: String): List<HlaSequence> {
        return readSequenceFiles({ x -> "${resourcesDir}/${x}_nuc.txt" }, { x -> x })
    }

    private fun readProteinFiles(resourcesDir: String): List<HlaSequence> {
        return readSequenceFiles({ x -> "${resourcesDir}/${x}_prot.txt" }, { x -> x.specificProteins() })
    }

    private fun readSequenceFiles(filenameSupplier: (Char) -> String, transform: (List<HlaSequence>) -> List<HlaSequence>): List<HlaSequence> {

        val aFile = filenameSupplier('A')
        val bFile = filenameSupplier('B')
        val cFile = filenameSupplier('C')

        val aSequence = transform(HlaSequenceFile.readFile(aFile).inflate())
        val bSequence = transform(HlaSequenceFile.readFile(bFile).inflate())
        val cSequence = transform(HlaSequenceFile.readFile(cFile).inflate())

        val result = mutableListOf<HlaSequence>()
        result.addAll(aSequence)
        result.addAll(bSequence)
        result.addAll(cSequence)

        val maxLength = result.map { it.sequence.length }.max()!!

        return result
                .filter { it.sequence.isNotEmpty() }
                .map { it.pad(maxLength) }
    }


    override fun close() {
        executorService.shutdown()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }

    private fun filterCandidatesOnNucleotides(minEvidence: Int, candidates: Collection<HlaSequence>, fragments: List<Fragment>, vararg nucleotides: Int): List<HlaSequence> {

        val reads = fragments
                .filter { it.containsAllNucleotides(*nucleotides) }
                .map { it.nucleotides(*nucleotides) }
                .groupingBy { it }
                .eachCount()
                .filter { it.value >= minEvidence }
                .keys
                .map { it.toCharArray() }


        return candidates.filter { it.consistentWith(nucleotides, reads) }

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


}
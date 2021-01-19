package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.lilac.candidates.Candidates
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.hla.HlaAlleleCoverageFactory
import com.hartwig.hmftools.lilac.hla.HlaAlleleCoverageFactory.Companion.coverageString
import com.hartwig.hmftools.lilac.hla.HlaComplex
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideFragmentFactory
import com.hartwig.hmftools.lilac.read.Fragment
import com.hartwig.hmftools.lilac.read.SAMRecordReader
import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.deflate
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.inflate
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.specificProteins
import org.apache.logging.log4j.LogManager
import java.util.concurrent.Executors

fun main(args: Array<String>) {
    LilacApplication().use { x -> x.run() }
}


class LilacApplication : AutoCloseable, Runnable {
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

    val resourcesDir = "/Users/jon/hmf/analysis/hla/resources"
    val outputDir = "/Users/jon/hmf/analysis/hla/resources"
//    val bamFile = "/Users/jon/hmf/analysis/hla/GIABvsSELFv004R.hla.bam"
    val bamFile = "/Users/jon/hmf/analysis/hla/COLO829v001R.hla.bam"
//    val bamFile = "/Users/jon/hmf/analysis/hla/COLO829v002R.hla.bam"
//    val bamFile = "/Users/jon/hmf/analysis/hla/COLO829v003R.hla.bam"


    override fun run() {
        logger.info("Starting LILAC")

        val aProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348, 364, 365)
        val bProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348, 362)
        val cProteinExonBoundaries = setOf(24, 114, 206, 298, 338, 349, 365, 366)
        val allProteinExonBoundaries = (aProteinExonBoundaries + bProteinExonBoundaries + cProteinExonBoundaries)
        val allNucleotideExonBoundaries = allProteinExonBoundaries.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }

        logger.info("Reading nucleotide files")
        val nucleotideSequences = readNucleotideFiles(resourcesDir)

        logger.info("Reading protein files")
        val aminoAcidSequences = readProteinFiles(resourcesDir)

        logger.info("Querying records from $bamFile")
        val rawNucleotideFragments = readFromBam(bamFile)
        val nucleotideFragmentFactory = NucleotideFragmentFactory(minBaseCount, rawNucleotideFragments, aProteinExonBoundaries, bProteinExonBoundaries, cProteinExonBoundaries)


        // Candidates
        val candidateFactory = Candidates(minBaseCount, 30, nucleotideSequences, aminoAcidSequences)
        val aCandidates = candidateFactory.candidates("A", aProteinExonBoundaries, nucleotideFragmentFactory.typeANucleotides())
        val bCandidates = candidateFactory.candidates("B", bProteinExonBoundaries, nucleotideFragmentFactory.typeBNucleotides())
        val cCandidates = candidateFactory.candidates("C", cProteinExonBoundaries, nucleotideFragmentFactory.typeCNucleotides())


        // Coverage
        val nucleotideCounts = SequenceCount.nucleotides(minBaseCount, nucleotideFragmentFactory.allNucleotides())
        val nucleotideHeterozygousLoci = nucleotideCounts.heterozygousIndices() intersect allNucleotideExonBoundaries
        val aminoAcidFragments = nucleotideFragmentFactory.allNucleotides().map { it.toAminoAcidFragment() }
        val aminoAcidCounts = SequenceCount.aminoAcids(minBaseCount, aminoAcidFragments)
        aminoAcidCounts.writeVertically("/Users/jon/hmf/analysis/hla/aminoacids.count.txt")


        val candidates = aCandidates + bCandidates + cCandidates
        val candidateAlleles = candidates.map { it.allele }
        val candidateAlleleSpecificProteins = candidateAlleles.map { it.specificProtein() }

        val aminoAcidCandidates = aminoAcidSequences.filter { it.allele in candidateAlleles }
        val nucleotideCandidates = nucleotideSequences.filter { it.allele.specificProtein() in candidateAlleleSpecificProteins }

        val coverageFactory = HlaAlleleCoverageFactory(aminoAcidFragments,
                aminoAcidCounts.heterozygousIndices(), aminoAcidCandidates,
                nucleotideHeterozygousLoci, nucleotideCandidates)


        logger.info("Calculating overall coverage")
        val groupCoverage = coverageFactory.groupCoverage(candidateAlleles)
        val confirmedGroups = groupCoverage.filter { it.uniqueCoverage > 0 }.sortedDescending()
        logger.info("... found ${confirmedGroups.size} uniquely identifiable groups: " + confirmedGroups.joinToString(", "))


        val proteinCoverage = coverageFactory.proteinCoverage(candidateAlleles)
        val confirmedProtein = proteinCoverage.filter { it.uniqueCoverage > 0 }.sortedDescending()
        logger.info("... found ${confirmedProtein.size} uniquely identifiable proteins: " + confirmedProtein.joinToString(", "))


        val complexes = HlaComplex.complexes(
                confirmedGroups.take(6).map { it.allele },
                confirmedProtein.take(6).map { it.allele },
                candidates.map { it.allele })

        logger.info("Calcuating coverage of ${complexes.size} complexes")

        println("TotalCoverage\tUniqueCoverage\tSharedCoverage\tAllele1\tAllele2\tAllele3\tAllele4\tAllele5\tAllele6")
        for (complex in complexes) {
            val complexCoverage = coverageFactory.proteinCoverage(complex.alleles)
            println(complexCoverage.coverageString())
        }


        val sequences = (candidates union (aminoAcidSequences.filter { it.allele == HlaAllele("C*03:04:01:01") })).map { HlaSequence(it.contig, it.sequence) }
        HlaSequenceFile.writeFile("/Users/jon/hmf/analysis/hla/candidates.inflate.txt", sequences)
        HlaSequenceFile.wipeFile("/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.writeBoundary(aProteinExonBoundaries, "/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.writeBoundary(bProteinExonBoundaries, "/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.writeBoundary(cProteinExonBoundaries, "/Users/jon/hmf/analysis/hla/candidates.deflate.txt")
        HlaSequenceFile.appendFile("/Users/jon/hmf/analysis/hla/candidates.deflate.txt", sequences.deflate())


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


    private fun readFromBam(bamFile: String): List<NucleotideFragment> {
        val transcripts = listOf(transcripts[HLA_A]!!, transcripts[HLA_B]!!, transcripts[HLA_C]!!)
        val reader = SAMRecordReader(1000, transcripts)
        val reads =  reader.readFromBam(bamFile)
        return NucleotideFragment.fromReads(minBaseQual, reads).filter { it.isNotEmpty() }
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

    private fun filterCandidatesOnNucleotides(minEvidence: Int, candidates: Collection<HlaSequence>, fragments: List<NucleotideFragment>, vararg nucleotides: Int): List<HlaSequence> {

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
package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.lilac.LilacApplication.Companion.logger
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.amino.AminoAcidFragmentPipeline
import com.hartwig.hmftools.lilac.candidates.Candidates
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceFactory
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceValidation
import com.hartwig.hmftools.lilac.hla.*
import com.hartwig.hmftools.lilac.hla.HlaComplexCoverage.Companion.writeToFile
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideGeneEnrichment
import com.hartwig.hmftools.lilac.read.SAMRecordReader
import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.inflate
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.specificProteins
import org.apache.commons.cli.*
import org.apache.logging.log4j.LogManager
import java.io.IOException
import java.util.concurrent.Executors
import kotlin.math.min

fun main(args: Array<String>) {

    @Throws(ParseException::class)
    fun createCommandLine(args: Array<String>, options: Options): CommandLine {
        val parser: CommandLineParser = DefaultParser()
        return parser.parse(options, args)
    }

    val options = LilacConfig.createOptions()
    try {
        val cmd = createCommandLine(args, options)
        val config = LilacConfig.createConfig(cmd)
        LilacApplication(config).use { x -> x.run() }
    } catch (e: IOException) {
        logger.warn(e)
    } catch (e: ParseException) {
        logger.warn(e)
        val formatter = HelpFormatter()
        formatter.printHelp("lilac", options)
    }

}


class LilacApplication(private val config: LilacConfig) : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
        const val HLA_A = "HLA-A"
        const val HLA_B = "HLA-B"
        const val HLA_C = "HLA-C"
    }

    private val startTime = System.currentTimeMillis()
    private val transcripts = HmfGenePanelSupplier.allGenesMap37()

    private val sample = config.sample
    private val minBaseQual = config.minBaseQual
    private val minEvidence = config.minEvidence
    private val minFragmentsPerAllele = config.minFragmentsPerAllele
    private val minFragmentsToRemoveSingle = config.minFragmentsToRemoveSingle
    private val minConfirmedUniqueCoverage = config.minConfirmedUniqueCoverage
    private val resourcesDir = config.resourceDir
    private val outputDir = config.outputDir
    private val bamFile = config.inputBam

    private val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("LILAC-%d").build()
    private val executorService = Executors.newFixedThreadPool(config.threads, namedThreadFactory)

    override fun run() {
        logger.info("Starting LILAC with parameters:")
        logger.info("     minFragmentsPerAllele = $minFragmentsPerAllele")
        logger.info("minFragmentsToRemoveSingle = $minFragmentsToRemoveSingle")
        logger.info("               minBaseQual = $minBaseQual")
        logger.info("               minEvidence = $minEvidence")
        logger.info("         minUniqueCoverage = $minConfirmedUniqueCoverage")

        val aProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348, 364, 365)
        val bProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348, 362)
        val cProteinExonBoundaries = setOf(24, 114, 206, 298, 338, 349, 365, 366)
        val allProteinExonBoundaries = (aProteinExonBoundaries + bProteinExonBoundaries + cProteinExonBoundaries)
        val allNucleotideExonBoundaries = allProteinExonBoundaries.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }

        // Context
        val hlaContextFactory = HlaContextFactory(aProteinExonBoundaries, bProteinExonBoundaries, cProteinExonBoundaries)
        val hlaAContext = hlaContextFactory.hlaA()
        val hlaBContext = hlaContextFactory.hlaB()
        val hlaCContext = hlaContextFactory.hlaC()

        logger.info("Querying records from $bamFile")
        val nucleotideGeneEnrichment = NucleotideGeneEnrichment(aProteinExonBoundaries, bProteinExonBoundaries, cProteinExonBoundaries)
        val rawNucleotideFragments = readFromBam(bamFile)
        val geneEnrichedNucleotideFragments = nucleotideGeneEnrichment.enrich(rawNucleotideFragments)
        val aminoAcidPipeline = AminoAcidFragmentPipeline(minBaseQual, minEvidence, geneEnrichedNucleotideFragments)
        val aFragments = aminoAcidPipeline.type(hlaAContext)
        val bFragments = aminoAcidPipeline.type(hlaBContext)
        val cFragments = aminoAcidPipeline.type(hlaCContext)

        logger.info("Reading nucleotide files")
        val nucleotideSequences = readNucleotideFiles(resourcesDir)

        logger.info("Reading protein files")
        val aminoAcidSequences = readProteinFiles(resourcesDir)

        // Phasing
        logger.info("Phasing bam records")
        val phasedEvidenceFactory = PhasedEvidenceFactory(minFragmentsPerAllele, config.minFragmentsToRemoveSingle, config.minEvidence)
        val aPhasedEvidence = phasedEvidenceFactory.evidence(hlaAContext, aFragments)
        val bPhasedEvidence = phasedEvidenceFactory.evidence(hlaBContext, bFragments)
        val cPhasedEvidence = phasedEvidenceFactory.evidence(hlaCContext, cFragments)

        // Validate phasing against expected sequences
        val expectedSequences = aminoAcidSequences.filter {  it.allele in config.expectedAlleles }
        PhasedEvidenceValidation.validateExpected("A", aPhasedEvidence, expectedSequences )
        PhasedEvidenceValidation.validateExpected("B", bPhasedEvidence, expectedSequences)
        PhasedEvidenceValidation.validateExpected("C", cPhasedEvidence,expectedSequences)

        // Candidates
        val candidateFactory = Candidates(config, nucleotideSequences, aminoAcidSequences)
        val aCandidates = candidateFactory.candidates(hlaAContext, aFragments, aPhasedEvidence)
        val bCandidates = candidateFactory.candidates(hlaBContext, bFragments, bPhasedEvidence)
        val cCandidates = candidateFactory.candidates(hlaCContext, cFragments, cPhasedEvidence)
        val candidates = aCandidates + bCandidates + cCandidates

        // Coverage
        val aminoAcidFragments = aminoAcidPipeline.combined(allProteinExonBoundaries)
        val nucleotideCounts = SequenceCount.nucleotides(minEvidence, aminoAcidFragments)
        val aminoAcidCounts = SequenceCount.aminoAcids(minEvidence, aminoAcidFragments)

        val nucleotideHeterozygousLoci = nucleotideCounts.heterozygousLoci() intersect allNucleotideExonBoundaries
        aminoAcidCounts.writeVertically("$outputDir/$sample.aminoacids.count.txt")
        nucleotideCounts.writeVertically("$outputDir/$sample.nucleotides.count.txt")


        val candidateAlleles = candidates.map { it.allele }
        val candidateAlleleSpecificProteins = candidateAlleles.map { it.specificProtein() }

        val aminoAcidCandidates = aminoAcidSequences.filter { it.allele in candidateAlleles }
        val nucleotideCandidates = nucleotideSequences.filter { it.allele.specificProtein() in candidateAlleleSpecificProteins }

        val coverageFactory = HlaComplexCoverageFactory(
                executorService,
                aminoAcidFragments,
                aminoAcidCounts.heterozygousLoci(),
                aminoAcidCandidates,
                nucleotideHeterozygousLoci,
                nucleotideCandidates)

        logger.info("Calculating overall allele[total,unique,shared,wide] coverage")
        val groupCoverage = coverageFactory.groupCoverage(candidateAlleles)
        val confirmedGroups = groupCoverage.alleles.filter { it.uniqueCoverage >= minConfirmedUniqueCoverage }.sortedDescending()
        val discardedGroups = groupCoverage.alleles.filter { it.uniqueCoverage in 1 until minConfirmedUniqueCoverage }.sortedDescending()
        logger.info("... found ${confirmedGroups.size} uniquely identifiable groups: " + confirmedGroups.joinToString(", "))
        if (discardedGroups.isNotEmpty()) {
            logger.info("... discarded ${discardedGroups.size} uniquely identifiable groups: " + discardedGroups.joinToString(", "))
        }

        val proteinCoverage = coverageFactory.proteinCoverage(candidateAlleles)
        val confirmedProtein = proteinCoverage.alleles.filter { it.uniqueCoverage >= minConfirmedUniqueCoverage }.sortedDescending()
        val discardedProtein = proteinCoverage.alleles.filter { it.uniqueCoverage in 1 until minConfirmedUniqueCoverage }.sortedDescending()
        logger.info("... found ${confirmedProtein.size} uniquely identifiable proteins: " + confirmedProtein.joinToString(", "))
        if (discardedProtein.isNotEmpty()) {
            logger.info("... discarded ${discardedProtein.size} uniquely identifiable proteins: " + discardedProtein.joinToString(", "))
        }

        val boundariesList = listOf(aProteinExonBoundaries, bProteinExonBoundaries, cProteinExonBoundaries)

        HlaSequenceFile.writeFile("$outputDir/$sample.candidates.inflate.txt", candidates)
        val template = aminoAcidSequences.filter { it.allele == HlaAllele("A*01:01:01:01") }.first()
        HlaSequenceFile.writeDeflatedFile("$outputDir/$sample.candidates.deflate.txt", boundariesList, template, candidates)

        val complexes = HlaComplex.complexes(
                confirmedGroups.take(6).map { it.allele },
                confirmedProtein.take(6).map { it.allele },
                candidates.map { it.allele })

        logger.info("Calculating coverage of ${complexes.size} complexes")
        val complexCoverage = coverageFactory.complexCoverage(complexes)
        if (complexCoverage.isNotEmpty()) {
            val topCoverage = complexCoverage[0]
            val minCoverage = min(topCoverage.totalCoverage * 0.99, topCoverage.totalCoverage - 10)
            val topComplexes = complexCoverage.filter { it.totalCoverage >= minCoverage }
            topComplexes.writeToFile("$outputDir/$sample.coverage.txt")

            logger.info(HlaComplexCoverage.header())
            for (topComplex in topComplexes) {
                logger.info(topComplex)
            }

            val winningAlleles = topCoverage.alleles.map { it.allele }
            val winningSequences = candidates.filter { candidate -> candidate.allele in winningAlleles }
            PhasedEvidenceValidation.validateAgainstFinalCandidates("A", aPhasedEvidence, winningSequences)
            PhasedEvidenceValidation.validateAgainstFinalCandidates("B", bPhasedEvidence, winningSequences)
            PhasedEvidenceValidation.validateAgainstFinalCandidates("C", cPhasedEvidence, winningSequences)
        }
    }

    private fun readFromBam(bamFile: String): List<NucleotideFragment> {
        val transcripts = listOf(transcripts[HLA_A]!!, transcripts[HLA_B]!!, transcripts[HLA_C]!!)
        val reader = SAMRecordReader(1000, config.refGenome, transcripts)
        val reads = reader.readFromBam(bamFile)
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

    private fun filterCandidatesOnExonBoundaryNucleotide(aminoAcidIndex: Int, minEvidence: Int, candidates: Collection<HlaSequence>, aminoAcidFragments: List<AminoAcidFragment>): List<HlaSequence> {
        val firstBaseCandidates = filterCandidatesOnNucleotides(minEvidence, candidates, aminoAcidFragments, aminoAcidIndex * 3)
        return filterCandidatesOnNucleotides(minEvidence, firstBaseCandidates, aminoAcidFragments, aminoAcidIndex * 3 + 1, aminoAcidIndex * 3 + 2)
    }

    private fun filterCandidatesOnExonBoundaries(exonBoundaries: Collection<Int>, minEvidence: Int, candidates: Collection<HlaSequence>, aminoAcidFragments: List<AminoAcidFragment>): List<HlaSequence> {
        var result = candidates.toList()
        for (exonBoundary in exonBoundaries) {
            result = filterCandidatesOnExonBoundaryNucleotide(exonBoundary, minEvidence, result, aminoAcidFragments)
        }

        return result
    }


}
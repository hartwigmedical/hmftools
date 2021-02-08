package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.lilac.LilacApplication.Companion.logger
import com.hartwig.hmftools.lilac.amino.AminoAcidFragmentPipeline
import com.hartwig.hmftools.lilac.candidates.Candidates
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceFactory
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceValidation
import com.hartwig.hmftools.lilac.hla.*
import com.hartwig.hmftools.lilac.hla.HlaComplexCoverage.Companion.writeToFile
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideFragmentFactory
import com.hartwig.hmftools.lilac.nuc.NucleotideGeneEnrichment
import com.hartwig.hmftools.lilac.read.SAMRecordReader
import com.hartwig.hmftools.lilac.seq.HlaSequence
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.inflate
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToFourDigit
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToSixDigit
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import com.hartwig.hmftools.lilac.seq.HlaSequenceLociFile
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


        val aProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348, 364)
        val bProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348)
        val cProteinExonBoundaries = setOf(24, 114, 206, 298, 338, 349, 365)
        val allProteinExonBoundaries = (aProteinExonBoundaries + bProteinExonBoundaries + cProteinExonBoundaries)
        val allNucleotideExonBoundaries = allProteinExonBoundaries.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }

        // Context
        val hlaContextFactory = HlaContextFactory(aProteinExonBoundaries, bProteinExonBoundaries, cProteinExonBoundaries)
        val hlaAContext = hlaContextFactory.hlaA()
        val hlaBContext = hlaContextFactory.hlaB()
        val hlaCContext = hlaContextFactory.hlaC()

        logger.info("Reading nucleotide files")
        val aNucleotideSequences = nucleotideLoci("${resourcesDir}/A_nuc.txt")
        val bNucleotideSequences = nucleotideLoci("${resourcesDir}/B_nuc.txt")
        val cNucleotideSequences = nucleotideLoci("${resourcesDir}/C_nuc.txt")
        val nucleotideSequences = aNucleotideSequences + bNucleotideSequences + cNucleotideSequences

        logger.info("Reading protein files")
        val aAminoAcidSequences = aminoAcidLoci("${resourcesDir}/A_prot.txt")
        val bAminoAcidSequences = aminoAcidLoci("${resourcesDir}/B_prot.txt")
        val cAminoAcidSequences = aminoAcidLoci("${resourcesDir}/C_prot.txt")
        val aminoAcidSequences = aAminoAcidSequences + bAminoAcidSequences + cAminoAcidSequences
        val aminoAcidSequencesWithInserts = aminoAcidSequences.filter { it.containsInserts() }
        val aminoAcidSequencesWithDeletes = aminoAcidSequences.filter { it.containsDeletes() }

        logger.info("Querying records from $bamFile")
        val nucleotideFragmentFactory = NucleotideFragmentFactory(config.minBaseQual, aminoAcidSequencesWithInserts, aminoAcidSequencesWithDeletes)
        val nucleotideGeneEnrichment = NucleotideGeneEnrichment(aProteinExonBoundaries, bProteinExonBoundaries, cProteinExonBoundaries)
        val rawNucleotideFragments = readFromBam(bamFile, nucleotideFragmentFactory)
        val geneEnrichedNucleotideFragments = nucleotideGeneEnrichment.enrich(rawNucleotideFragments)
        val aminoAcidPipeline = AminoAcidFragmentPipeline(minBaseQual, minEvidence, geneEnrichedNucleotideFragments)
        val aFragments = aminoAcidPipeline.type(hlaAContext)
        val bFragments = aminoAcidPipeline.type(hlaBContext)
        val cFragments = aminoAcidPipeline.type(hlaCContext)


        // Phasing
        logger.info("Phasing bam records")
        val phasedEvidenceFactory = PhasedEvidenceFactory(minFragmentsPerAllele, config.minFragmentsToRemoveSingle, config.minEvidence)
        val aPhasedEvidence = phasedEvidenceFactory.evidence(hlaAContext, aFragments)
        val bPhasedEvidence = phasedEvidenceFactory.evidence(hlaBContext, bFragments)
        val cPhasedEvidence = phasedEvidenceFactory.evidence(hlaCContext, cFragments)

        // Validate phasing against expected sequences
        val expectedSequences = aminoAcidSequences.filter { it.allele in config.expectedAlleles }
        PhasedEvidenceValidation.validateExpected("A", aPhasedEvidence, expectedSequences)
        PhasedEvidenceValidation.validateExpected("B", bPhasedEvidence, expectedSequences)
        PhasedEvidenceValidation.validateExpected("C", cPhasedEvidence, expectedSequences)

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
        val candidateAlleleSpecificProteins = candidateAlleles.map { it.asFourDigit() }

        val aminoAcidCandidates = aminoAcidSequences.filter { it.allele in candidateAlleles }
        val nucleotideCandidates = nucleotideSequences.filter { it.allele.asFourDigit() in candidateAlleleSpecificProteins }

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

        val candidatesAfterConfirmedGroups = candidateAlleles.filterWithConfirmedGroups(confirmedGroups.map { it.allele })
        val proteinCoverage = coverageFactory.proteinCoverage(candidatesAfterConfirmedGroups)
        val confirmedProtein = proteinCoverage.alleles.filter { it.uniqueCoverage >= minConfirmedUniqueCoverage }.sortedDescending()
        val discardedProtein = proteinCoverage.alleles.filter { it.uniqueCoverage in 1 until minConfirmedUniqueCoverage }.sortedDescending()
        logger.info("... found ${confirmedProtein.size} uniquely identifiable proteins: " + confirmedProtein.joinToString(", "))
        if (discardedProtein.isNotEmpty()) {
            logger.info("... discarded ${discardedProtein.size} uniquely identifiable proteins: " + discardedProtein.joinToString(", "))
        }

        val template = aminoAcidSequences.filter { it.allele == HlaAllele("A*01:01:01:01") }.first()
        val writtenCandidates = (candidates + expectedSequences + template).distinct().sortedBy { it.allele }
        HlaSequenceLociFile.write("$outputDir/$sample.candidates.deflate.txt", writtenCandidates)

        val complexes = HlaComplex.complexes(
                confirmedGroups.take(6).map { it.allele },
                confirmedProtein.take(6).map { it.allele },
                candidates.map { it.allele })

        logger.info("Calculating coverage of ${complexes.size} complexes")

        if (expectedSequences.isNotEmpty()) {
            val expectedCoverage = coverageFactory.proteinCoverage(expectedSequences.map { it.allele })
            logger.info("Expected coverage: $expectedCoverage")
        }

        val complexCoverage = coverageFactory.complexCoverage(complexes)
        if (complexCoverage.isNotEmpty()) {
            val topCoverage = complexCoverage[0]
            val minCoverage = min(topCoverage.totalCoverage * 0.99, topCoverage.totalCoverage - 10.0)
            val topComplexes = complexCoverage.filter { it.totalCoverage >= minCoverage }
            topComplexes.writeToFile("$outputDir/$sample.coverage.txt")

            logger.info(HlaComplexCoverage.header())
            for (topComplex in topComplexes) {
                logger.info(topComplex)
            }

            val complexesWithSameTopScore = complexCoverage.filter { it.totalCoverage == topCoverage.totalCoverage }.count()
            logger.info("${config.sample} - $complexesWithSameTopScore WINNERS, WINNING ALLELES: ${topCoverage.alleles.map { it.allele }}")

            val winningAlleles = topCoverage.alleles.map { it.allele }
            val winningSequences = candidates.filter { candidate -> candidate.allele in winningAlleles }
            PhasedEvidenceValidation.validateAgainstFinalCandidates("A", aPhasedEvidence, winningSequences)
            PhasedEvidenceValidation.validateAgainstFinalCandidates("B", bPhasedEvidence, winningSequences)
            PhasedEvidenceValidation.validateAgainstFinalCandidates("C", cPhasedEvidence, winningSequences)
        }
    }

    private fun readFromBam(bamFile: String, nucleotideFragmentFactory: NucleotideFragmentFactory): List<NucleotideFragment> {
        val transcripts = listOf(transcripts[HLA_A]!!, transcripts[HLA_B]!!, transcripts[HLA_C]!!)
        val reader = SAMRecordReader(1000, config.refGenome, transcripts, nucleotideFragmentFactory)
        return reader.readFromBam(bamFile)
    }

    fun unwrapFile(fileName: String) {
        val input = HlaSequenceFile.readFile(fileName)
        val output = input.reduceToFourDigit()
        val inflated = output.map { it.inflate(input[0].rawSequence) }
        val deflated = inflated.map { it.deflate(inflated[0].rawSequence) }
        HlaSequenceFile.writeFile(fileName.replace(".txt", ".unwrapped.txt"), output)
        HlaSequenceFile.writeFile(fileName.replace(".txt", ".deflated.txt"), deflated)
    }

    private fun readNucleotideFiles(resourcesDir: String): List<HlaSequence> {
        return readSequenceFiles({ x -> "${resourcesDir}/${x}_nuc.txt" }, { x -> x })
    }

    private fun readProteinFiles(resourcesDir: String): List<HlaSequence> {
        return readSequenceFiles({ x -> "${resourcesDir}/${x}_prot.txt" }, { x -> x.reduceToFourDigit() })
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


    private fun nucleotideLoci(inputFilename: String): List<HlaSequenceLoci> {
        val sequences = HlaSequenceFile.readFile(inputFilename)
//                .filter { it.allele.gene in setOf("A", "B") ||  it.allele in listOf(HlaAllele("C*01:02:01:01"), HlaAllele("C*17:01:01:02")) }
                .reduceToSixDigit()
        return HlaSequenceLoci.create(sequences)
                .filter { it.sequences.isNotEmpty() }
    }

    private fun aminoAcidLoci(inputFilename: String): List<HlaSequenceLoci> {
        val sequences = HlaSequenceFile.readFile(inputFilename)
                .reduceToFourDigit()
                .map { if (it.rawSequence.endsWith("X")) it else it.copyWithAdditionalSequence("X") }

        return HlaSequenceLoci.create(sequences)
                .filter { it.sequences.isNotEmpty() }
    }


    override fun close() {
        executorService.shutdown()
        logger.info("Finished in ${(System.currentTimeMillis() - startTime) / 1000} seconds")
    }

    private fun List<HlaAllele>.filterWithConfirmedGroups(confirmedGroups: List<HlaAllele>): List<HlaAllele> {
        val a = confirmedGroups.filter { it.gene == "A" }
        val b = confirmedGroups.filter { it.gene == "B" }
        val c = confirmedGroups.filter { it.gene == "C" }
        val map = mapOf(Pair("A", a), Pair("B", b), Pair("C", c))

        return this.filter { map[it.gene]!!.size < 2 || map[it.gene]!!.contains(it.asAlleleGroup()) }

    }

}
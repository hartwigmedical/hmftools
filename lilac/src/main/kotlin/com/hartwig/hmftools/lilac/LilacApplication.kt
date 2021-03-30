package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.lilac.LilacApplication.Companion.logger
import com.hartwig.hmftools.lilac.amino.AminoAcidFragmentPipeline
import com.hartwig.hmftools.lilac.candidates.Candidates
import com.hartwig.hmftools.lilac.cna.HlaCopyNumber
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage
import com.hartwig.hmftools.lilac.coverage.HlaComplex
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverage.Companion.writeToFile
import com.hartwig.hmftools.lilac.coverage.HlaComplexCoverageFactory
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceFactory
import com.hartwig.hmftools.lilac.evidence.PhasedEvidenceValidation
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.hla.HlaContextFactory
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideFragmentFactory
import com.hartwig.hmftools.lilac.nuc.NucleotideGeneEnrichment
import com.hartwig.hmftools.lilac.out.HlaOut
import com.hartwig.hmftools.lilac.qc.*
import com.hartwig.hmftools.lilac.read.FragmentAlleles
import com.hartwig.hmftools.lilac.read.SAMRecordReader
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToFourDigit
import com.hartwig.hmftools.lilac.seq.HlaSequenceFile.reduceToSixDigit
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import com.hartwig.hmftools.lilac.seq.HlaSequenceLociFile
import org.apache.commons.cli.*
import org.apache.logging.log4j.LogManager
import java.io.IOException
import java.util.concurrent.Executors
import kotlin.system.exitProcess

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
        exitProcess(1)
    } catch (e: ParseException) {
        logger.warn(e)
        val formatter = HelpFormatter()
        formatter.printHelp("lilac", options)
        exitProcess(1)
    }

}


class LilacApplication(private val config: LilacConfig) : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
        const val HLA_A = "HLA-A"
        const val HLA_B = "HLA-B"
        const val HLA_C = "HLA-C"

        val DEFLATE_TEMPLATE = HlaAllele("A*01:01")
        val EXCLUDED_ALLELES = setOf(HlaAllele("A*01:81"), HlaAllele("A*01:237"))

        val A_EXON_BOUNDARIES = setOf(24, 114, 206, 298, 337, 348, 364)
        val B_EXON_BOUNDARIES = setOf(24, 114, 206, 298, 337, 348)
        val C_EXON_BOUNDARIES = setOf(24, 114, 206, 298, 338, 349, 365)

        val ALL_PROTEIN_EXON_BOUNDARIES = (A_EXON_BOUNDARIES + B_EXON_BOUNDARIES + C_EXON_BOUNDARIES)
        val ALL_NUCLEOTIDE_EXON_BOUNDARIES = ALL_PROTEIN_EXON_BOUNDARIES.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }
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

    private val namedThreadFactory = ThreadFactoryBuilder().setNameFormat("LILAC-%d").build()
    private val executorService = Executors.newFixedThreadPool(config.threads, namedThreadFactory)
    private val nucleotideGeneEnrichment = NucleotideGeneEnrichment(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES)

    override fun run() {
        logger.info("Starting LILAC with parameters:")
        logger.info("    minBaseQual = $minBaseQual")
        logger.info("    minEvidence = $minEvidence")
        logger.info("    minUniqueCoverage = $minConfirmedUniqueCoverage")
        logger.info("    minFragmentsPerAllele = $minFragmentsPerAllele")
        logger.info("    minFragmentsToRemoveSingle = $minFragmentsToRemoveSingle")

        // Context
        val hlaContextFactory = HlaContextFactory(A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES)
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

        logger.info("Querying records from reference bam ${config.referenceBam}")
        val nucleotideFragmentFactory = NucleotideFragmentFactory(config.minBaseQual, aminoAcidSequencesWithInserts, aminoAcidSequencesWithDeletes)
        val transcripts = listOf(transcripts[HLA_A]!!, transcripts[HLA_B]!!, transcripts[HLA_C]!!)
        val bamReader = SAMRecordReader(1000, config.refGenome, transcripts, nucleotideFragmentFactory)
        val referenceNucleotideFragments = bamReader.readFromBam(config.referenceBam).enrichGenes()

        val tumorNucleotideFragments = if (config.tumorBam.isNotEmpty()) {
            logger.info("Querying records from tumor bam ${config.tumorBam}")
            bamReader.readFromBam(config.tumorBam).enrichGenes()
        } else {
            listOf()
        }

        logger.info("Enriching reference bam records")
        val aminoAcidPipeline = AminoAcidFragmentPipeline(config, referenceNucleotideFragments, tumorNucleotideFragments)
        val aCandidateFragments = aminoAcidPipeline.referenceCandidateFragments(hlaAContext)
        val bCandidateFragments = aminoAcidPipeline.referenceCandidateFragments(hlaBContext)
        val cCandidateFragments = aminoAcidPipeline.referenceCandidateFragments(hlaCContext)

        // Phasing
        logger.info("Phasing bam records")
        val phasedEvidenceFactory = PhasedEvidenceFactory(minFragmentsPerAllele, config.minFragmentsToRemoveSingle, config.minEvidence)
        val aPhasedEvidence = phasedEvidenceFactory.evidence(hlaAContext, aCandidateFragments)
        val bPhasedEvidence = phasedEvidenceFactory.evidence(hlaBContext, bCandidateFragments)
        val cPhasedEvidence = phasedEvidenceFactory.evidence(hlaCContext, cCandidateFragments)

        // Validate phasing against expected sequences
        val expectedSequences = aminoAcidSequences.filter { it.allele.asFourDigit() in config.expectedAlleles }
        PhasedEvidenceValidation.validateExpected("A", aPhasedEvidence, expectedSequences)
        PhasedEvidenceValidation.validateExpected("B", bPhasedEvidence, expectedSequences)
        PhasedEvidenceValidation.validateExpected("C", cPhasedEvidence, expectedSequences)

        // Candidates
        val candidateFactory = Candidates(config, nucleotideSequences, aminoAcidSequences)
        val aCandidates = candidateFactory.candidates(hlaAContext, aCandidateFragments, aPhasedEvidence)
        val bCandidates = candidateFactory.candidates(hlaBContext, bCandidateFragments, bPhasedEvidence)
        val cCandidates = candidateFactory.candidates(hlaCContext, cCandidateFragments, cPhasedEvidence)
        val candidates = aCandidates + bCandidates + cCandidates

        // Coverage
        val referenceCoverageFragments = aminoAcidPipeline.referenceCoverageFragments()
        val referenceAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, referenceCoverageFragments)
        val referenceAminoAcidHeterozygousLoci = referenceAminoAcidCounts.heterozygousLoci()
        val referenceNucleotideCounts = SequenceCount.nucleotides(minEvidence, referenceCoverageFragments)
        val referenceNucleotideHeterozygousLoci = referenceNucleotideCounts.heterozygousLoci() intersect ALL_NUCLEOTIDE_EXON_BOUNDARIES

        val candidateAlleles = candidates.map { it.allele }
        val candidateAlleleSpecificProteins = candidateAlleles.map { it.asFourDigit() }
        val candidateAminoAcidSequences = aminoAcidSequences.filter { it.allele in candidateAlleles }
        val candidateNucleotideSequences = nucleotideSequences.filter { it.allele.asFourDigit() in candidateAlleleSpecificProteins }

        val referenceFragmentAlleles = FragmentAlleles.create(referenceCoverageFragments,
                referenceAminoAcidHeterozygousLoci, candidateAminoAcidSequences, referenceNucleotideHeterozygousLoci, candidateNucleotideSequences)
        val coverageFactory = HlaComplexCoverageFactory(config, executorService)


        logger.info("Identifying uniquely identifiable groups and proteins [total,unique,shared,wide]")
        val groupCoverage = HlaComplexCoverageFactory.groupCoverage(referenceFragmentAlleles, candidateAlleles)
        val confirmedGroups = groupCoverage.confirmUnique()
        val discardedGroups = groupCoverage.alleleCoverage.filter { it.uniqueCoverage > 0 && it !in confirmedGroups }.sortedDescending()
        if (confirmedGroups.isNotEmpty()) {
            logger.info("    confirmed ${confirmedGroups.size} unique groups: " + confirmedGroups.joinToString(", "))
        } else {
            logger.info("    confirmed 0 unique groups")
        }
        if (discardedGroups.isNotEmpty()) {
            logger.info("    found ${discardedGroups.size} insufficiently unique groups: " + discardedGroups.joinToString(", "))
        }

        val candidatesAfterConfirmedGroups = candidateAlleles.filterWithConfirmedGroups(confirmedGroups.map { it.allele })
        val proteinCoverage = HlaComplexCoverageFactory.proteinCoverage(referenceFragmentAlleles, candidatesAfterConfirmedGroups)
        val confirmedProtein = proteinCoverage.confirmUnique()
        val discardedProtein = proteinCoverage.alleleCoverage.filter { it.uniqueCoverage > 0 && it !in confirmedProtein }.sortedDescending()
        if (confirmedProtein.isNotEmpty()) {
            logger.info("    confirmed ${confirmedProtein.size} unique proteins: " + confirmedProtein.joinToString(", "))
        } else {
            logger.info("    confirmed 0 unique proteins")
        }
        if (discardedProtein.isNotEmpty()) {
            logger.info("    found ${discardedProtein.size} insufficiently unique proteins: " + discardedProtein.joinToString(", "))
        }

        val complexes = HlaComplex.complexes(confirmedGroups.alleles(), confirmedProtein.alleles(), candidates.map { it.allele })
        logger.info("Calculating coverage of ${complexes.size} complexes")

        val referenceRankedComplexes = coverageFactory.rankedComplexCoverage(referenceFragmentAlleles, complexes)
        if (referenceRankedComplexes.isEmpty()) {
            logger.fatal("Failed to calculate complex coverage")
            exitProcess(1)
        }

        if (expectedSequences.isNotEmpty()) {
            val expectedCoverage = HlaComplexCoverageFactory.proteinCoverage(referenceFragmentAlleles, expectedSequences.map { it.allele })
            logger.info("Expected allele coverage: $expectedCoverage")
        }

        val winningReferenceCoverage = referenceRankedComplexes[0].expandToSixAlleles()
        val winningAlleles = winningReferenceCoverage.alleleCoverage.alleles()
        val winningSequences = candidates.filter { candidate -> candidate.allele in winningAlleles }

        logger.info(HlaComplexCoverage.header())
        for (rankedComplex in referenceRankedComplexes) {
            logger.info(rankedComplex)
        }

        logger.info("${config.sample} - REF - ${referenceRankedComplexes.size} CANDIDATES, WINNING ALLELES: ${winningReferenceCoverage.alleleCoverage.map { it.allele }}")

        val aminoAcidQC = AminoAcidQC.create(winningSequences, referenceAminoAcidCounts)
        val haplotypeQC = HaplotypeQC.create(3, winningSequences, aPhasedEvidence + bPhasedEvidence + cPhasedEvidence, referenceAminoAcidCounts)
        val bamQC = BamQC.create(bamReader)
        val coverageQC = CoverageQC.create(referenceFragmentAlleles.size, winningReferenceCoverage)
        val lilacQC = LilacQC(aminoAcidQC, bamQC, coverageQC, haplotypeQC)

        logger.info("QC Stats:")
        logger.info("    ${lilacQC.header().joinToString(",")}")
        logger.info("    ${lilacQC.body().joinToString(",")}")

        lilacQC.writefile("$outputDir/$sample.qc.txt")


        val winningTumorCoverage: HlaComplexCoverage
        val winningTumorCopyNumber: HlaCopyNumber

        if (config.tumorBam.isNotEmpty()) {
            logger.info("Calculating tumor coverage of winning alleles")

            val tumorFragmentAlleles = FragmentAlleles.create(aminoAcidPipeline.tumorCoverageFragments(),
                    referenceAminoAcidHeterozygousLoci, candidateAminoAcidSequences, referenceNucleotideHeterozygousLoci, candidateNucleotideSequences)
            winningTumorCoverage = HlaComplexCoverageFactory.proteinCoverage(tumorFragmentAlleles, winningAlleles).expandToSixAlleles()

            logger.info("Calculating tumor copy number of winning alleles")
            winningTumorCopyNumber = HlaCopyNumber.alleleCopyNumber(config.geneCopyNumberFile, winningTumorCoverage)
        } else {
            winningTumorCoverage = HlaComplexCoverage.create(listOf())
            winningTumorCopyNumber = HlaCopyNumber.alleleCopyNumber(config.geneCopyNumberFile, winningTumorCoverage)
        }

        val output = HlaOut(winningReferenceCoverage, winningTumorCoverage, winningTumorCopyNumber)
        output.write("$outputDir/$sample.lilac.txt")

        logger.info("Writing output to $outputDir")
        val deflatedSequenceTemplate = aminoAcidSequences.first { it.allele == DEFLATE_TEMPLATE }
        val candidateToWrite = (candidates + expectedSequences + deflatedSequenceTemplate).distinct().sortedBy { it.allele }
        HlaSequenceLociFile.write("$outputDir/$sample.candidates.sequences.txt", A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES, candidateToWrite)
        referenceRankedComplexes.writeToFile("$outputDir/$sample.candidates.coverage.txt")
        referenceAminoAcidCounts.writeVertically("$outputDir/$sample.aminoacids.count.txt")
        referenceNucleotideCounts.writeVertically("$outputDir/$sample.nucleotides.count.txt")

    }


    private fun nucleotideLoci(inputFilename: String): List<HlaSequenceLoci> {
        val sequences = HlaSequenceFile.readFile(inputFilename)
                .filter { it.allele !in EXCLUDED_ALLELES }
                .reduceToSixDigit()
        return HlaSequenceLoci.create(sequences)
                .filter { it.sequences.isNotEmpty() }
    }

    private fun aminoAcidLoci(inputFilename: String): List<HlaSequenceLoci> {
        val sequences = HlaSequenceFile.readFile(inputFilename)
                .filter { it.allele !in EXCLUDED_ALLELES }
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

    private fun HlaComplexCoverage.confirmUnique(): List<HlaAlleleCoverage> {
        val unique = this.alleleCoverage.filter { it.uniqueCoverage >= config.minConfirmedUniqueCoverage }.sortedDescending()
        val a = unique.filter { it.allele.gene == "A" }.take(2)
        val b = unique.filter { it.allele.gene == "B" }.take(2)
        val c = unique.filter { it.allele.gene == "C" }.take(2)

        return (a + b + c).sortedDescending()
    }


    private fun List<HlaAlleleCoverage>.alleles(): List<HlaAllele> {
        return this.map { it.allele }
    }

    private fun List<NucleotideFragment>.enrichGenes(): List<NucleotideFragment> {
        return nucleotideGeneEnrichment.enrich(this)
    }

}
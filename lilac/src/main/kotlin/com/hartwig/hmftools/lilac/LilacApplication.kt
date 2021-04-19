package com.hartwig.hmftools.lilac

import com.google.common.util.concurrent.ThreadFactoryBuilder
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier
import com.hartwig.hmftools.common.hla.HlaFiles
import com.hartwig.hmftools.common.utils.version.VersionInfo
import com.hartwig.hmftools.common.variant.VariantContextDecorator
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
import com.hartwig.hmftools.lilac.variant.LilacVCF
import com.hartwig.hmftools.lilac.variant.SomaticAlleleCoverage
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount.Companion.addVariant
import com.hartwig.hmftools.lilac.variant.SomaticVariants
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
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
        LilacApplication(cmd, config).use { x -> x.run() }
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


class LilacApplication(private val cmd: CommandLine, private val config: LilacConfig) : AutoCloseable, Runnable {
    companion object {
        val logger = LogManager.getLogger(this::class.java)
        const val HLA_A = "HLA-A"
        const val HLA_B = "HLA-B"
        const val HLA_C = "HLA-C"
        val HLA_GENES = setOf(HLA_A, HLA_B, HLA_C)
        val VERSION = VersionInfo("lilac.version")

        val DEFLATE_TEMPLATE = HlaAllele("A*01:01")
        val EXCLUDED_ALLELES = setOf(
                HlaAllele("A*01:81"),
                HlaAllele("A*01:237"),
                HlaAllele("A*33:191"),
                HlaAllele("A*11:353"),
                HlaAllele("A*30:95"),
                HlaAllele("A*30:136"),
                HlaAllele("A*31:135"))

        val A_EXON_BOUNDARIES = setOf(24, 114, 206, 298, 337, 348, 364)
        val B_EXON_BOUNDARIES = setOf(24, 114, 206, 298, 337, 348)
        val C_EXON_BOUNDARIES = setOf(24, 114, 206, 298, 338, 349, 365)

        val ALL_PROTEIN_EXON_BOUNDARIES = (A_EXON_BOUNDARIES + B_EXON_BOUNDARIES + C_EXON_BOUNDARIES)
        val ALL_NUCLEOTIDE_EXON_BOUNDARIES = ALL_PROTEIN_EXON_BOUNDARIES.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }

        private val transcripts = HmfGenePanelSupplier.allGenesMap37()
        val HLA_TRANSCRIPTS = listOf(
                transcripts[HLA_A]!!,
                transcripts[HLA_B]!!,
                transcripts[HLA_C]!!
        )
        val LOCI_POSITION = LociPosition(HLA_TRANSCRIPTS)
    }

    private val startTime = System.currentTimeMillis()
    private val sample = config.sample
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
        logger.info("    sample = ${config.sample}")
        logger.info("    minBaseQual = ${config.minBaseQual}")
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
        val nucleotideFragmentFactory = NucleotideFragmentFactory(config.minBaseQual, aminoAcidSequencesWithInserts, aminoAcidSequencesWithDeletes, LOCI_POSITION)
        val tumorBamReader = SAMRecordReader(config.tumorBam, config.refGenome, HLA_TRANSCRIPTS, nucleotideFragmentFactory)
        val referenceBamReader = SAMRecordReader(config.referenceBam, config.refGenome, HLA_TRANSCRIPTS, nucleotideFragmentFactory)
        val referenceNucleotideFragments = referenceBamReader.readFromBam().enrichGenes()

        val tumorNucleotideFragments = if (config.tumorBam.isNotEmpty()) {
            logger.info("Querying records from tumor bam ${config.tumorBam}")
            tumorBamReader.readFromBam().enrichGenes()
        } else {
            listOf()
        }

        // Enrich bam records
        val aminoAcidPipeline = AminoAcidFragmentPipeline(config, referenceNucleotideFragments, tumorNucleotideFragments)
        val aCandidateFragments = aminoAcidPipeline.referencePhasingFragments(hlaAContext)
        val bCandidateFragments = aminoAcidPipeline.referencePhasingFragments(hlaBContext)
        val cCandidateFragments = aminoAcidPipeline.referencePhasingFragments(hlaCContext)

        // Un-phased Candidates
        val candidateFactory = Candidates(config, nucleotideSequences, aminoAcidSequences)
        val aUnphasedCandidates = candidateFactory.unphasedCandidates(hlaAContext, aCandidateFragments)
        val bUnphasedCandidates = candidateFactory.unphasedCandidates(hlaBContext, bCandidateFragments)
        val cUnphasedCandidates = candidateFactory.unphasedCandidates(hlaCContext, cCandidateFragments)
        val allUnphasedCandidates = (aUnphasedCandidates + bUnphasedCandidates + cUnphasedCandidates)

        // Phasing
        val phasedEvidenceFactory = PhasedEvidenceFactory(config)
        val aPhasedEvidence = phasedEvidenceFactory.evidence(hlaAContext, aCandidateFragments)
        val bPhasedEvidence = phasedEvidenceFactory.evidence(hlaBContext, bCandidateFragments)
        val cPhasedEvidence = phasedEvidenceFactory.evidence(hlaCContext, cCandidateFragments)

        // Validate phasing against expected sequences
        val expectedSequences = aminoAcidSequences.filter { it.allele.asFourDigit() in config.expectedAlleles }
        PhasedEvidenceValidation.validateExpected("A", aPhasedEvidence, expectedSequences)
        PhasedEvidenceValidation.validateExpected("B", bPhasedEvidence, expectedSequences)
        PhasedEvidenceValidation.validateExpected("C", cPhasedEvidence, expectedSequences)

        // Phased Candidates
        val aCandidates = candidateFactory.phasedCandidates(hlaAContext, aUnphasedCandidates, aPhasedEvidence)
        val bCandidates = candidateFactory.phasedCandidates(hlaBContext, bUnphasedCandidates, bPhasedEvidence)
        val cCandidates = candidateFactory.phasedCandidates(hlaCContext, cUnphasedCandidates, cPhasedEvidence)
        val phasedCandidateAlleles = (aCandidates + bCandidates + cCandidates)

        val missingStopLossAlleles = config.stopLossRecoveryAlleles.filter { it !in phasedCandidateAlleles }
        val stopLossRecovery = if (nucleotideFragmentFactory.containsStopLossHlaC() && missingStopLossAlleles.isNotEmpty()) {
            logger.info("Recovering stop loss candidate: " + missingStopLossAlleles.joinToString(","))
            missingStopLossAlleles
        } else {
            listOf()
        }

        // Common Candidate Recovery
        logger.info("Recovering common un-phased candidates:")
        val recoveredAlleles = config.commonAlleles
                .filter { it !in phasedCandidateAlleles && it in allUnphasedCandidates } + stopLossRecovery

        val candidateAlleles = phasedCandidateAlleles + recoveredAlleles
        val candidateSequences = aminoAcidSequences.filter { it.allele in candidateAlleles }
        if (recoveredAlleles.isNotEmpty()) {
            logger.info("    recovered ${recoveredAlleles.size} common candidate alleles: " + recoveredAlleles.joinToString(", "))
        } else {
            logger.info("    recovered 0 common candidate alleles")
        }

        // Coverage
        val referenceCoverageFragments = aminoAcidPipeline.referenceCoverageFragments()
        val referenceAminoAcidCounts = SequenceCount.aminoAcids(minEvidence, referenceCoverageFragments)
        val referenceAminoAcidHeterozygousLoci = referenceAminoAcidCounts.heterozygousLoci()
        val referenceNucleotideCounts = SequenceCount.nucleotides(minEvidence, referenceCoverageFragments)
        val referenceNucleotideHeterozygousLoci = referenceNucleotideCounts.heterozygousLoci() intersect ALL_NUCLEOTIDE_EXON_BOUNDARIES

        val candidateAlleleSpecificProteins = candidateAlleles.map { it.asFourDigit() }
        val candidateAminoAcidSequences = aminoAcidSequences.filter { it.allele in candidateAlleles }
        val candidateNucleotideSequences = nucleotideSequences.filter { it.allele.asFourDigit() in candidateAlleleSpecificProteins }

        val referenceFragmentAlleles = FragmentAlleles.create(referenceCoverageFragments,
                referenceAminoAcidHeterozygousLoci, candidateAminoAcidSequences, referenceNucleotideHeterozygousLoci, candidateNucleotideSequences)

        // Complexes
        val complexes = HlaComplex.complexes(config, referenceFragmentAlleles, candidateAlleles, recoveredAlleles)

        logger.info("Calculating coverage of ${complexes.size} complexes")
        val coverageFactory = HlaComplexCoverageFactory(config)
        val referenceRankedComplexes = coverageFactory.rankedComplexCoverage(executorService, referenceFragmentAlleles, complexes, recoveredAlleles)
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
        val winningSequences = candidateSequences.filter { candidate -> candidate.allele in winningAlleles }.toSet()

        logger.info(HlaComplexCoverage.header())
        for (rankedComplex in referenceRankedComplexes) {
            logger.info(rankedComplex)
        }

        logger.info("${config.sample} - REF - ${referenceRankedComplexes.size} CANDIDATES, WINNING ALLELES: ${winningReferenceCoverage.alleleCoverage.map { it.allele }}")

        val winningTumorCoverage: HlaComplexCoverage
        val winningTumorCopyNumber: List<HlaCopyNumber>
        val somaticVariants: List<VariantContextDecorator>
        var somaticCodingCount = SomaticCodingCount.create(winningAlleles)

        if (config.tumorBam.isNotEmpty()) {
            logger.info("Calculating tumor coverage of winning alleles")
            val tumorFragmentAlleles = FragmentAlleles.create(aminoAcidPipeline.tumorCoverageFragments(),
                    referenceAminoAcidHeterozygousLoci, candidateAminoAcidSequences, referenceNucleotideHeterozygousLoci, candidateNucleotideSequences)
            winningTumorCoverage = HlaComplexCoverageFactory.proteinCoverage(tumorFragmentAlleles, winningAlleles).expandToSixAlleles()

            logger.info("Calculating tumor copy number of winning alleles")
            winningTumorCopyNumber = HlaCopyNumber.alleleCopyNumber(winningAlleles, config.geneCopyNumberFile, winningTumorCoverage)

            // SOMATIC VARIANTS
            somaticVariants = SomaticVariants(config).readSomaticVariants()
            if (somaticVariants.isNotEmpty()) {
                logger.info("Calculating somatic variant allele coverage")
                val lilacVCF = LilacVCF("${config.outputFilePrefix}.lilac.somatic.vcf.gz", config.somaticVcf).writeHeader()

                val somaticCoverageFactory = SomaticAlleleCoverage(config, referenceAminoAcidHeterozygousLoci, LOCI_POSITION, somaticVariants, winningSequences)
                for (variant in somaticVariants) {
                    val variantCoverage = somaticCoverageFactory.alleleCoverage(variant, tumorBamReader)
                    val variantAlleles = variantCoverage.alleles().toSet()
                    logger.info("    $variant -> $variantCoverage")
                    lilacVCF.writeVariant(variant.context(), variantAlleles)
                    somaticCodingCount = somaticCodingCount.addVariant(variant, variantAlleles)
                }
            }
        } else {
            winningTumorCoverage = HlaComplexCoverage.create(listOf())
            winningTumorCopyNumber = HlaCopyNumber.alleleCopyNumber(winningAlleles)
            somaticVariants = listOf()
        }

        val output = HlaOut.create(winningReferenceCoverage, winningTumorCoverage, winningTumorCopyNumber, somaticCodingCount)

        // QC
        logger.info("Calculating QC Statistics")
        val somaticVariantQC = SomaticVariantQC.create(somaticVariants.size, somaticCodingCount)
        val aminoAcidQC = AminoAcidQC.create(winningSequences, referenceAminoAcidCounts)
        val haplotypeQC = HaplotypeQC.create(3, winningSequences, aPhasedEvidence + bPhasedEvidence + cPhasedEvidence, referenceAminoAcidCounts)
        val bamQC = BamQC.create(referenceBamReader)
        val coverageQC = CoverageQC.create(referenceNucleotideFragments.size, winningReferenceCoverage)
        val lilacQC = LilacQC.create(aminoAcidQC, bamQC, coverageQC, haplotypeQC, somaticVariantQC)

        logger.info("QC Stats:")
        logger.info("    ${lilacQC.header().joinToString(",")}")
        logger.info("    ${lilacQC.body().joinToString(",")}")

        logger.info("Writing output to $outputDir")
        val outputFile = "${config.outputFilePrefix}.lilac.txt"
        val outputQCFile = "${config.outputFilePrefix}.lilac.qc.txt"

        output.write(outputFile)
        lilacQC.writefile(outputQCFile)
        val deflatedSequenceTemplate = aminoAcidSequences.first { it.allele == DEFLATE_TEMPLATE }
        val candidateToWrite = (candidateSequences + expectedSequences + deflatedSequenceTemplate).distinct().sortedBy { it.allele }
        HlaSequenceLociFile.write("$outputDir/$sample.candidates.sequences.txt", A_EXON_BOUNDARIES, B_EXON_BOUNDARIES, C_EXON_BOUNDARIES, candidateToWrite)
        referenceRankedComplexes.writeToFile("$outputDir/$sample.candidates.coverage.txt")
        referenceAminoAcidCounts.writeVertically("$outputDir/$sample.candidates.aminoacids.txt")
        referenceNucleotideCounts.writeVertically("$outputDir/$sample.candidates.nucleotides.txt")

        if (DatabaseAccess.hasDatabaseConfig(cmd)) {
            logger.info("Writing output to DB")
            val dbAccess = DatabaseAccess.databaseAccess(cmd, true)
            val type = HlaFiles.type(outputFile, outputQCFile);
            val typeDetails = HlaFiles.typeDetails(outputFile)
            dbAccess.writeHla(sample, type, typeDetails)
        }
    }

    private fun nucleotideLoci(inputFilename: String): List<HlaSequenceLoci> {
        val sequences = HlaSequenceFile.readFile(inputFilename)
                .filter { it.allele !in EXCLUDED_ALLELES }
                .reduceToSixDigit()
        return HlaSequenceLoci.create(sequences)
                .filter { it.sequences.isNotEmpty() }
    }


    private fun createSynonymous(template: HlaSequenceLoci): HlaSequenceLoci {
        val modSequences = template.sequences.toMutableList()
        modSequences[1011] = "A"
        modSequences[1012] = "G"
        modSequences[1013] = "T"

        return HlaSequenceLoci(HlaAllele("C*04:82:01"), modSequences)
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


    private fun List<HlaAlleleCoverage>.alleles(): List<HlaAllele> {
        return this.map { it.allele }
    }

    private fun List<NucleotideFragment>.enrichGenes(): List<NucleotideFragment> {
        return nucleotideGeneEnrichment.enrich(this)
    }

}
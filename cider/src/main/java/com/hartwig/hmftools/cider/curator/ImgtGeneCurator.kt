package com.hartwig.hmftools.cider.curator

import com.beust.jcommander.JCommander
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.UnixStyleUsageFormatter
import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.cider.CiderConstants.BLAST_REF_GENOME_VERSION
import com.hartwig.hmftools.cider.CiderConstants.BLASTN_PRIMARY_ASSEMBLY_NAME
import com.hartwig.hmftools.cider.annotation.ImgtSequenceFile
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.ANCHOR_MISMATCH_MAX
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.BLASTN_EVALUE_CUTOFF
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.BLASTN_MAX_MISMATCH
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IGKDEL_IMGT_SEQ
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IGKINTR_IMGT_SEQ
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IMGT_ANCHOR_LENGTH
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IMGT_V_ANCHOR_INDEX
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.FASTA_REF_CONTEXT
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.REF_CONTEXT_CHECK
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.REF_CONTEXT_CHECK_MISMATCH_MAX
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.SPECIES
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.jAnchorSignatures
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.liftOverBlacklist
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.cider.genes.IgTcrFunctionality
import com.hartwig.hmftools.cider.genes.IgTcrGene
import com.hartwig.hmftools.cider.genes.IgTcrGene.Companion.toCommonIgTcrGene
import com.hartwig.hmftools.cider.genes.IgTcrRegion
import com.hartwig.hmftools.cider.genes.VJGeneType
import com.hartwig.hmftools.common.bam.CigarUtils.getPositionFromReadIndex
import com.hartwig.hmftools.common.blastn.BlastnMatch
import com.hartwig.hmftools.common.cider.IgTcrGeneFile
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache
import com.hartwig.hmftools.common.gene.GeneData
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.utils.Doubles
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import htsjdk.samtools.CigarElement
import htsjdk.samtools.liftover.LiftOver
import htsjdk.samtools.util.Interval
import htsjdk.samtools.util.SequenceUtil.reverseComplement
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import java.io.File
import java.nio.file.Paths
import java.util.Comparator
import kotlin.math.min
import kotlin.system.exitProcess

const val ANCHOR_DNA_LENGTH: Int = 30

data class ProcessedGeneAllele(
    val imgt: ImgtGeneAllele,
    val locationV38: GenomicLocation?,
    val locationV37: GenomicLocation?,
    val anchorSequence: String?,
    val anchorLocationV38: GenomicLocation?,
    val anchorLocationV37: GenomicLocation?,
)

data class LocationInfo(
    val location: GenomicLocation,
    val cigar: List<CigarElement>?,
)

data class AnchorInfo(
    val sequence: String,
    val location: GenomicLocation?
)

// This utility processes IMGT IG/TCR genes into format CIDER can use.
// It uses BLASTN combined with ensembl to find the genomic location of each IG/TCR genes.
// The result is then lifted over to V37.
// Download the IMGT sequences from the IMGT website.
//
// example usage:
// java -cp VjTemplateGeneWriter
// -imgt IMGT_genes.fasta
// -blast /tools/ncbi-blast
// -blast_db /data/blast_db
// -ensembl_data_dir /data/resources/public/ensembl_data_cache/38
// -liftover_chain /data/hg38ToHg19.over.chain.gz
// -ref_genome_v38 /data/resources/bucket/reference_genome/38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
// -ref_genome_v37 /data/resources/bucket/reference_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta
// -output_dir ./output
//
// Format for the fasta line is F+ORF+in-frame P nucleotide sequences with IMGT gaps
// See IMGT doc: https://www.imgt.org/genedb/doc
//
// The liftover chain file can be downloaded from
// https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
class ImgtGeneCurator
{
    @Parameter(names = ["-imgt"], required = true, description = "Input IMGT fasta file")
    lateinit var inputImgtFasta: String

    @Parameter(names = ["-blast"], required = true, description = "Location of blast installation")
    lateinit var blast: String

    @Parameter(names = ["-blast_db"], required = true, description = "Location of blast database")
    lateinit var blastDb: String

    @Parameter(names = ["-" + EnsemblDataCache.ENSEMBL_DATA_DIR], required = true, description = EnsemblDataCache.ENSEMBL_DATA_DIR_CFG)
    lateinit var ensemblDataDir: String

    @Parameter(names = ["-liftover_chain"], required = true, description = "38 to 37 liftover chain file")
    lateinit var liftOverChainFile: String

    @Parameter(names = ["-ref_genome_v38"], required = true, description = "V38 Reference genome fasta file")
    lateinit var refGenomeV38: String

    @Parameter(names = ["-ref_genome_v37"], required = true, description = "V37 Reference genome fasta file")
    lateinit var refGenomeV37: String

    @Parameter(names = ["-output_dir"], required = true, description = "Output directory")
    lateinit var outputDir: String

    @Parameter(names = ["-threads"], description = "Number of threads")
    var threadCount = 1

    @Parameter(names = ["-workdir"], description = "Number of threads")
    lateinit var workdir: String

    lateinit var refGenomeSourceV38: RefGenomeSource
    lateinit var refGenomeSourceV37: RefGenomeSource
    lateinit var ensemblDataCache: EnsemblDataCache
    lateinit var liftOverHtsjdk: LiftOver
    lateinit var liftOverHmf: GenomeLiftoverCache

    fun run(): Int
    {
        refGenomeSourceV38 = loadRefGenome(refGenomeV38)
        refGenomeSourceV37 = loadRefGenome(refGenomeV37)

        ensemblDataCache = EnsemblDataCache(ensemblDataDir, BLAST_REF_GENOME_VERSION)
        val ensemblLoadOk = ensemblDataCache.load(true)
        if (!ensemblLoadOk)
        {
            sLogger.error("Ensembl data cache load failed")
            throw RuntimeException("Ensembl data cache load failed")
        }

        liftOverHtsjdk = LiftOver(File(liftOverChainFile))
        liftOverHmf = GenomeLiftoverCache(true)

        val geneAlleles = loadGeneAlleles(inputImgtFasta)
        val processedGeneAlleles = processAlleles(geneAlleles)
        val (ciderGenesV38, ciderGenesV37) = createCiderGenes(processedGeneAlleles)

        sLogger.info("validating V38 anchor locations")
        validateAnchorLocations(ciderGenesV38, refGenomeSourceV38)
        sLogger.info("validating V37 anchor locations")
        validateAnchorLocations(ciderGenesV37, refGenomeSourceV37)

        val outputGenesV38 = Paths.get(outputDir, "igtcr_gene.38.tsv").toString()
        IgTcrGeneFile.write(outputGenesV38, ciderGenesV38.map { o -> toCommonIgTcrGene(o) })
        val outputGenesV37 = Paths.get(outputDir, "igtcr_gene.37.tsv").toString()
        IgTcrGeneFile.write(outputGenesV37, ciderGenesV37.map { o -> toCommonIgTcrGene(o) })

        sLogger.info("Writing V38 allele fasta")
        val outputFastaV38 = Paths.get(outputDir, "igtcr_gene.38.fasta").toString()
        writeVDJFasta(outputFastaV38, processedGeneAlleles, true)
        sLogger.info("Writing V37 allele fasta")
        val outputFastaV37 = Paths.get(outputDir, "igtcr_gene.37.fasta").toString()
        writeVDJFasta(outputFastaV37, processedGeneAlleles, false)

        return 0
    }

    private fun processAlleles(alleles: List<ImgtGeneAllele>): List<ProcessedGeneAllele>
    {
        val alleleLocationsV38 = findAlleleLocations(alleles)
        return alleles.mapIndexed { index, allele -> processAllele(allele, alleleLocationsV38[index]) }
    }

    private fun processAllele(allele: ImgtGeneAllele, alleleLocationInfoV38: LocationInfo?): ProcessedGeneAllele
    {
        val region = allele.region!!

        val anchorInfo = when (region)
        {
            IgTcrRegion.V_REGION -> findVAnchor(allele, alleleLocationInfoV38)
            IgTcrRegion.J_REGION -> findJAnchor(allele, alleleLocationInfoV38)
            else -> null
        }

        val (alleleLocationV37, anchorLocationV37) = convertAlleleLocationsTo37(allele, anchorInfo?.sequence, alleleLocationInfoV38?.location, anchorInfo?.location)

        return ProcessedGeneAllele(allele, alleleLocationInfoV38?.location, alleleLocationV37, anchorInfo?.sequence, anchorInfo?.location, anchorLocationV37)
    }

    private fun findAlleleLocations(alleles: List<ImgtGeneAllele>): List<LocationInfo?>
    {
        // V / D / J gene use blast, constant genes use ensembl
        // reason is that ensembl is easier, and we do not need the alt locations for the constant genes. Another reason
        // is that we do not need to be very precise with the location of the anchor for constant genes
        // if we use ensembl for V / J gene, we need to use the fasta file to validate the anchor location is precise

        val locations = MutableList<LocationInfo?>(alleles.size) { null }

        val nonConstantAlleles = alleles.withIndex().filter { (_, allele) -> allele.region != IgTcrRegion.CONSTANT }
        val nonConstantAlleleLocations = blastForAlleleLocation(nonConstantAlleles.map { it.value })
        nonConstantAlleleLocations.withIndex().forEach { (indexInNonconstant, location) ->
            locations[nonConstantAlleles[indexInNonconstant].index] = location
        }

        val constantAlleles = alleles.withIndex().filter { (_, gene) -> gene.region == IgTcrRegion.CONSTANT }
        for ((index, allele) in constantAlleles)
        {
            val ensemblGene = ensemblDataCache.getGeneDataByName(allele.geneName) ?: continue
            locations[index] = LocationInfo(toGenomicLocation(ensemblGene), null)
        }

        for ((index, allele) in alleles.withIndex())
        {
            ImgtGeneCuratorSettings.getGenomicLocationOverrides(allele.geneName)
                ?.let { override ->
                    sLogger.info("gene: {}, using location override: {}", allele.geneAllele, override)
                    locations[index] = override
                }
        }

        return locations
    }

    private fun blastForAlleleLocation(alleles: List<ImgtGeneAllele>): List<LocationInfo?>
    {
        val blastnResults = runBlastn(
            "imgt", blast, blastDb,
            alleles.map { it -> it.sequenceWithoutGaps },
            workdir, threadCount, BLASTN_EVALUE_CUTOFF)

        val locations = ArrayList<LocationInfo?>()
        for ((alleleIndex, allele) in alleles.withIndex())
        {
            // Filter by full match and not too many mismatches, then
            // We sort the matches by
            // 1. match quality
            // 2. primary assembly
            val matches = blastnResults[alleleIndex]
                .filter { m ->
                    m.numMismatch <= BLASTN_MAX_MISMATCH &&
                            m.querySeqLen - (m.queryAlignEnd - m.queryAlignStart + 1) <= BLASTN_MAX_MISMATCH }
                .sortedWith(Comparator.comparingDouble { m: BlastnMatch -> m.expectedValue }
                    .thenComparingInt { m: BlastnMatch -> if (m.subjectTitle.endsWith(BLASTN_PRIMARY_ASSEMBLY_NAME)) 0 else 1 })

            // find the gene in the ensembl
            val ensemblGene = ensemblDataCache.getGeneDataByName(allele.geneName)
            var bestMatch: BlastnMatch? = null

            for (match in matches)
            {
                val matchLocation = blastnMatchtoGenomicLocation(match)

                if (matchLocation == null)
                {
                    continue
                }

                if (bestMatch == null)
                {
                    bestMatch = match
                }

                else if (Doubles.equal(bestMatch.bitScore, match.bitScore))
                {
                    // if they have same score, we want to choose the one that overlaps with ensembl gene
                    if (ensemblGene != null &&
                        ensemblGene.Chromosome == matchLocation.chromosome &&
                        ensemblGene.forwardStrand() == (matchLocation.strand == Strand.FORWARD) &&
                        ensemblGene.GeneStart <= matchLocation.posEnd &&
                        ensemblGene.GeneEnd >= matchLocation.posStart)
                    {
                        bestMatch = match
                        sLogger.debug("gene: {}, using match that overlaps with ensembl", allele.geneAllele)
                    }
                }
            }

            var location: LocationInfo? = null
            if (bestMatch == null)
            {
                if (ensemblGene != null)
                {
                    sLogger.error("gene: {}, no full match yet has ensembl", allele.geneAllele)

                    if (allele.region == IgTcrRegion.D_REGION)
                    {
                        // use ensembl for D region
                        location = LocationInfo(toGenomicLocation(ensemblGene), null)
                    }
                }
                sLogger.info("gene: {}, no full match", allele.geneAllele)
            }
            else
            {
                val genomicLocation = matchToQueryGenomicLocation(bestMatch)
                location = LocationInfo(genomicLocation, bestMatch.cigar)
                sLogger.info(
                    "gene: {}, gene loc: {}, contig={}, start={}, end={}, cigar: {}",
                    allele.geneAllele, genomicLocation,
                    bestMatch.subjectTitle, bestMatch.subjectAlignStart, bestMatch.subjectAlignEnd, bestMatch.cigar)
            }
            locations.add(location)
        }
        return locations
    }

    private fun convertAlleleLocationsTo37(
        allele: ImgtGeneAllele, anchorSequence: String?, alleleLocationV38: GenomicLocation?, anchorLocationV38: GenomicLocation?
    ): Pair<GenomicLocation?, GenomicLocation?>
    {
        var alleleLocationV37 = alleleLocationV38?.let(this::convertGenomicLocationTo37)

        var anchorLocationV37 = anchorLocationV38?.let(this::convertGenomicLocationTo37)

        // apply blacklist
        if (allele.geneName in liftOverBlacklist)
        {
            sLogger.info("gene: {} in liftover blacklist, clearing v37 genomic location", allele.geneName)
            alleleLocationV37 = null
            anchorLocationV37 = null
        }

        // some anchor lengths are not the same, but we cannot do much about them
        if (anchorLocationV37 != null && anchorLocationV37.baseLength() != anchorLocationV38!!.baseLength())
        {
            // decide whether to change anchor start or anchor end
            // the logic is that for V gene we need to preserve the TGT at the end of anchor
            // for J gene we need to preserve the TGG at the start of anchor
            if ((allele.region == IgTcrRegion.V_REGION) == (anchorLocationV37.strand == Strand.FORWARD))
            {
                sLogger.warn("gene: {}, strand: {}, different base lengths between v38({}) and v37({}), anchor length: {}, changing anchor start",
                    allele.geneAllele, anchorLocationV37.strand, anchorLocationV38.baseLength(), anchorLocationV37.baseLength(),
                    anchorSequence!!.length)

                anchorLocationV37 = anchorLocationV37.copy(posStart = anchorLocationV37.posEnd - anchorSequence.length + 1)
            }
            else
            {
                sLogger.warn("gene: {}, strand: {}, different base lengths between v38({}) and v37({}), anchor length: {}, changing anchor end",
                    allele.geneAllele, anchorLocationV37.strand, anchorLocationV38.baseLength(), anchorLocationV37.baseLength(),
                    anchorSequence!!.length)

                anchorLocationV37 = anchorLocationV37.copy(posEnd = anchorLocationV37.posStart + anchorSequence.length - 1)
            }
        }

        return Pair(alleleLocationV37, anchorLocationV37)
    }

    // To convert from v38 to v37, we use two methods. First we use the Hartwig GenomicLiftOverCache, and then
    // we use the htsjdk Liftover with downloaded Chain file.
    private fun convertGenomicLocationTo37(locationV38: GenomicLocation): GenomicLocation?
    {
        if (locationV38.inPrimaryAssembly)
        {
            // try the HMF one
            val genomicLocV37 = convertGenomicLocationTo37Hmf(locationV38)

            if (genomicLocV37 != null)
            {
                return genomicLocV37
            }
            else
            {
                val interval38 = Interval(locationV38.chromosome, locationV38.posStart, locationV38.posEnd, locationV38.strand == Strand.REVERSE, "")
                val interval37 = liftOverHtsjdk.liftOver(interval38)

                if (interval37 != null)
                {
                    sLogger.info("HMF cannot convert but htsjdk liftover can: v38({}) v37({})", locationV38, genomicLocV37)
                    return locationV38.copy(chromosome = RefGenomeVersion.V37.versionedChromosome(locationV38.chromosome),
                        posStart = interval37.start, posEnd = interval37.end, strand = if (interval37.isPositiveStrand) Strand.FORWARD else Strand.REVERSE)
                }
            }
        }
        return null
    }

    private fun convertGenomicLocationTo37Hmf(locationV38: GenomicLocation): GenomicLocation?
    {
        if (locationV38.inPrimaryAssembly)
        {
            val v37PosStart: Int = liftOverHmf.convertPositionTo37(locationV38.chromosome, locationV38.posStart)
            val v37PosEnd: Int = liftOverHmf.convertPositionTo37(locationV38.chromosome, locationV38.posEnd)

            if (v37PosStart != -1 && v37PosEnd != -1)
            {
                return if (v37PosEnd < v37PosStart)
                {
                    // different strand
                    locationV38.copy(chromosome = RefGenomeVersion.V37.versionedChromosome(locationV38.chromosome),
                        posStart = v37PosEnd, posEnd = v37PosStart, strand = locationV38.strand.opposite)
                }
                else
                {
                    locationV38.copy(chromosome = RefGenomeVersion.V37.versionedChromosome(locationV38.chromosome),
                        posStart = v37PosStart, posEnd = v37PosEnd)
                }
            }
        }
        return null
    }

    private fun createCiderGenes(processedAlleles: List<ProcessedGeneAllele>): Pair<List<IgTcrGene>, List<IgTcrGene>>
    {
        val genesV38 = processedAlleles.map { createCiderGene(it, true) }
        val genesV37 = processedAlleles.map { createCiderGene(it, false) }
        return Pair(genesV38, genesV37)
    }

    private fun createCiderGene(processedAllele: ProcessedGeneAllele, isV38: Boolean): IgTcrGene
    {
        return IgTcrGene(
            processedAllele.imgt.geneName,
            processedAllele.imgt.allele,
            processedAllele.imgt.region!!,
            processedAllele.imgt.functionality,
            if (isV38) processedAllele.locationV38 else processedAllele.locationV37,
            processedAllele.anchorSequence,
            if (isV38) processedAllele.anchorLocationV38 else processedAllele.anchorLocationV37,
        )
    }

    private fun writeVDJFasta(path: String, alleles: List<ProcessedGeneAllele>, isV38: Boolean)
    {
        File(path).printWriter().use { file ->
            alleles.forEach { allele ->
                if (allele.imgt.region?.isVDJ ?: false)
                {
                    val sequence = getAlleleSequenceWithContext(allele, isV38)
                    file.print(">${sequence.fastaLabel}\n${sequence.sequenceWithRef}\n")
                }
            }
        }
    }

    private fun getAlleleSequenceWithContext(allele: ProcessedGeneAllele, isV38: Boolean): ImgtSequenceFile.Sequence
    {
        val alleleSeq = allele.imgt.sequenceWithoutGaps
        val location = if (isV38) allele.locationV38 else allele.locationV37
        val isForward = location?.strand == Strand.FORWARD

        val seq: String
        val refBeforeLength: Int
        val refAfterLength: Int

        // If it's not in the primary assembly then the chromosome from Blastn won't index the ref genome.
        if (location == null || !location.inPrimaryAssembly)
        {
            seq = alleleSeq
            refBeforeLength = 0
            refAfterLength = 0
        }
        else
        {
            val alleleSeqAligned = if (isForward) alleleSeq else reverseComplement(alleleSeq)

            val refGenome = if (isV38) refGenomeSourceV38 else refGenomeSourceV37
            val chromosome = getChromosomeForRefGenome(location.chromosome, refGenome)
            // When we get the reference sequence surrounding the allele, check a few bases overlapping the allele to ensure it lines up.
            val ref1 = refGenome.getBaseString(
                chromosome,
                location.posStart - FASTA_REF_CONTEXT,
                location.posStart - 1 + REF_CONTEXT_CHECK
            ) ?: ""
            val ref2 = refGenome.getBaseString(
                chromosome,
                location.posEnd + 1 - REF_CONTEXT_CHECK,
                location.posEnd + FASTA_REF_CONTEXT
            ) ?: ""
            val refStart = ref1.substring(FASTA_REF_CONTEXT)
            val alleleStart = alleleSeqAligned.substring(0, REF_CONTEXT_CHECK)
            val refEnd = ref2.substring(0, REF_CONTEXT_CHECK)
            val alleleEnd = alleleSeqAligned.substring(alleleSeqAligned.length - REF_CONTEXT_CHECK)
            val startMismatches = countMismatches(alleleStart, refStart)
            val endMismatches = countMismatches(alleleEnd, refEnd)
            if (startMismatches > REF_CONTEXT_CHECK_MISMATCH_MAX || endMismatches > REF_CONTEXT_CHECK_MISMATCH_MAX)
            {
                sLogger.error(
                    "Gene allele ref mismatch: allele={} alleleStart={} refStart={} alleleEnd={} refEnd={}",
                    allele.imgt.geneAllele, alleleStart, refStart, alleleEnd, refEnd
                )
                seq = alleleSeq
                refBeforeLength = 0
                refAfterLength = 0
            }
            else
            {
                val refBefore = ref1.substring(0, ref1.length - REF_CONTEXT_CHECK)
                val refAfter = ref2.substring(REF_CONTEXT_CHECK)
                val seqAligned = refBefore + alleleSeqAligned + refAfter
                // Note the final sequence is always on the same strand as the original IMGT resource.
                seq = if (isForward) seqAligned else reverseComplement(seqAligned)
                refBeforeLength = if (isForward) refBefore.length else refAfter.length
                refAfterLength = if (isForward) refAfter.length else refBefore.length
            }
        }

        return ImgtSequenceFile.Sequence(allele.imgt.geneName, allele.imgt.allele, seq, refBeforeLength, refAfterLength)
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(ImgtGeneCurator::class.java)

        @JvmStatic
        fun main(args: Array<String>)
        {
            val app = ImgtGeneCurator()
            val commander = JCommander.newBuilder()
                .addObject(app)
                .build()

            // use unix style formatter
            commander.usageFormatter = UnixStyleUsageFormatter(commander)
            commander.parameterDescriptionComparator = DeclaredOrderParameterComparator(app.javaClass)

            try
            {
                commander.parse(*args)
                exitProcess(app.run())
            }
            catch (paramException: ParameterException)
            {
                println("${paramException.message}")
                commander.usage()
                exitProcess(1)
            }
        }

        private fun loadGeneAlleles(imgtFastaPath: String): List<ImgtGeneAllele>
        {
            val imgtAlleles = readGeneDataFromFasta(imgtFastaPath)

            // Add IGKINTR and IGKDEL.
            val alleles = imgtAlleles.toMutableList()
            alleles.add(ImgtGeneAllele(
                geneName = VJGeneType.IGKINTR, allele = "01", species = SPECIES, functionality = IgTcrFunctionality.ORF,
                region = IgTcrRegion.V_REGION, sequenceWithGaps = IGKINTR_IMGT_SEQ, partial = false
            ))
            alleles.add(ImgtGeneAllele(
                geneName = VJGeneType.IGKDEL, allele = "01", species = SPECIES, functionality = IgTcrFunctionality.ORF,
                region = IgTcrRegion.J_REGION, sequenceWithGaps = IGKDEL_IMGT_SEQ, partial = false
            ))

            val filteredAlleles = alleles.filter {
                it.species == SPECIES && it.region != null
            }

            filteredAlleles.forEach { sLogger.info("Gene allele: $it") }

            return filteredAlleles
        }

        private fun findVAnchor(allele: ImgtGeneAllele, alleleLocationInfo: LocationInfo?): AnchorInfo?
        {
            // some sequences are longer
            val seqWithGaps = allele.sequenceWithGaps

            if (seqWithGaps.length < IMGT_V_ANCHOR_INDEX)
            {
                sLogger.log(
                    if (allele.functionality == IgTcrFunctionality.FUNCTIONAL && !allele.partial) Level.ERROR else Level.INFO,
                    "Cannot find V anchor, sequence too short. {}", allele.geneAllele)
                return null
            }

            val anchor = seqWithGaps.substring(IMGT_V_ANCHOR_INDEX, min(IMGT_V_ANCHOR_INDEX + IMGT_ANCHOR_LENGTH, seqWithGaps.length))

            // if anchor is too short we remove
            if (anchor.length < min(IMGT_ANCHOR_LENGTH, ANCHOR_DNA_LENGTH))
            {
                // skip this one
                sLogger.log(
                    if (allele.functionality == IgTcrFunctionality.FUNCTIONAL && !allele.partial) Level.ERROR else Level.INFO,
                    "V anchor:{} too short for {}", anchor, allele.geneAllele
                )
                return null
            }

            // anchor cannot contain .
            if (anchor.contains("."))
            {
                sLogger.error("V anchor: {} contains \".\", gene: {}", anchor, allele.geneAllele)
                return null
            }

            val seqWithoutGaps = allele.sequenceWithoutGaps

            // now we find the anchor index again using the one without gap
            val anchorIndex = seqWithoutGaps.indexOf(anchor)

            var anchorLocation: GenomicLocation? = null

            val alleleLocation = alleleLocationInfo?.location
            if (alleleLocation != null && anchorIndex >= 0)
            {
                // Work based on the end of the anchor for V genes, because it's more important that lines up correctly in the case of indels.

                val anchorEndInSeq = anchorIndex + anchor.length - 1
                val anchorEndInAligned = if (alleleLocation.strand == Strand.FORWARD) anchorEndInSeq else seqWithoutGaps.length - anchorEndInSeq - 1
                val anchorEndOffsetFromStart = alleleLocationInfo.cigar
                    ?.let { getPositionFromReadIndex(0, it, anchorEndInAligned) }
                    ?: anchorEndInAligned

                anchorLocation = if (alleleLocation.strand == Strand.FORWARD)
                {
                    val anchorEndOffsetFromEnd = alleleLocation.baseLength() - anchorEndOffsetFromStart - 1
                    alleleLocation.copy(
                        posStart = alleleLocation.posEnd - anchorEndOffsetFromEnd - (anchor.length - 1),
                        posEnd = alleleLocation.posEnd - anchorEndOffsetFromEnd
                    )
                } else
                {
                    alleleLocation.copy(
                        posStart = alleleLocation.posStart + anchorEndOffsetFromStart,
                        posEnd = alleleLocation.posStart + anchorEndOffsetFromStart + (anchor.length - 1)
                    )
                }
            }

            val aaSeq = Codons.aminoAcidFromBases(anchor)

            // v gene
            sLogger.info(
                "V gene: {}, anchor: {}, offset from start: {}, anchor AA: {}",
                allele.geneAllele,
                anchor,
                anchorIndex,
                aaSeq
            )

            return AnchorInfo(anchor, anchorLocation)
        }

        private fun findJAnchor(allele: ImgtGeneAllele, alleleLocationInfo: LocationInfo?): AnchorInfo?
        {
            // J gene rules
            // 30 base sequence starting with TGGGG (W) or TTTG and TTCG (F)
            val seqWithGaps = allele.sequenceWithGaps
            val anchorIndex: Int = seqWithGaps.indexOfAny(jAnchorSignatures(allele.geneName))

            if (anchorIndex <= 0)
            {
                sLogger.log(
                    if (allele.functionality == IgTcrFunctionality.FUNCTIONAL && !allele.partial) Level.ERROR else Level.INFO,
                    "J gene: {} cannot find anchor", allele.geneAllele)
                return null
            }

            val anchor = seqWithGaps.substring(anchorIndex, min(anchorIndex + IMGT_ANCHOR_LENGTH, seqWithGaps.length))

            // cannot contain .
            if (anchor.contains("."))
            {
                sLogger.error("J gene: {}, anchor({}) contains .", allele.geneAllele, anchor)
                return null
            }

            var anchorLocation: GenomicLocation? = null

            val alleleLocation = alleleLocationInfo?.location
            if (alleleLocation != null)
            {
                // TODO: this should be indel aware
                // TODO: why is this using the index in the sequence with gaps?
                anchorLocation = if (alleleLocation.strand == Strand.FORWARD)
                {
                    alleleLocation.copy(
                        posStart = alleleLocation.posStart + anchorIndex,
                        posEnd = alleleLocation.posStart + anchorIndex + anchor.length - 1
                    )
                } else
                {
                    alleleLocation.copy(
                        posStart = alleleLocation.posEnd - anchorIndex - anchor.length + 1,
                        posEnd = alleleLocation.posEnd - anchorIndex
                    )
                }
            }

            val aaSeq = Codons.aminoAcidFromBases(anchor)

            sLogger.info("J gene: {}, anchor: {}, offset from start: {}, anchor AA: {}", allele.geneAllele, anchor, anchorIndex, aaSeq)

            return AnchorInfo(anchor, anchorLocation)
        }

        private fun toGenomicLocation(ensemblGene: GeneData) : GenomicLocation
        {
            return GenomicLocation(ensemblGene.Chromosome, ensemblGene.GeneStart,
                ensemblGene.GeneEnd, Strand.valueOf(ensemblGene.Strand.toInt()))
        }

        internal fun matchToQueryGenomicLocation(match: BlastnMatch) : GenomicLocation
        {
            // we need to correct for the ends to make sure things align properly
            val startExtend = match.queryAlignStart - 1
            val endExtend = match.querySeqLen - match.queryAlignEnd
            val matchGenomicLoc = blastnMatchtoGenomicLocation(match)!!

            if (startExtend == 0 && endExtend == 0)
            {
                return matchGenomicLoc
            }

            return if (matchGenomicLoc.strand == Strand.FORWARD)
            {
                matchGenomicLoc.copy(posStart = matchGenomicLoc.posStart - startExtend, posEnd = matchGenomicLoc.posEnd + endExtend)
            }
            else
            {
                matchGenomicLoc.copy(posStart = matchGenomicLoc.posStart - endExtend, posEnd = matchGenomicLoc.posEnd + startExtend)
            }
        }

        private fun validateAnchorLocations(igTcrGeneList: List<IgTcrGene>, refGenome: RefGenomeSource)
        {
            var failed = false
            for (gene in igTcrGeneList)
            {
                if (gene.anchorLocation != null && gene.anchorLocation.inPrimaryAssembly)
                {
                    if (!validateAgainstRefGenome(gene.anchorSequence!!, gene.anchorLocation, refGenome, ANCHOR_MISMATCH_MAX))
                    {
                        sLogger.error("gene: {} anchor location: {} does not match anchor seq: {}",
                            gene.geneAllele, gene.anchorLocation, gene.anchorSequence)
                        failed = true
                    }
                }
            }
            if (failed)
            {
                //throw RuntimeException("Invalid anchor locations")
            }
        }

        // validate sequence against the ref genome file to make sure we got it right
        fun validateAgainstRefGenome(seq: String, genomicLocation: GenomicLocation, refGenome: RefGenomeSource, maxMismatches: Int): Boolean
        {
            val refForwardSeq = refGenome.getBaseString(getChromosomeForRefGenome(genomicLocation.chromosome, refGenome), genomicLocation.posStart, genomicLocation.posEnd)
            val refGenomeSeq = if (genomicLocation.strand == Strand.FORWARD) refForwardSeq else reverseComplement(refForwardSeq)

            if (refGenomeSeq.length != seq.length)
            {
                sLogger.warn("validation failed: seq({}) and ref genome seq({} of {}) length mismatch", seq, refGenomeSeq, genomicLocation)
                return false
            }

            val mismatches = countMismatches(seq, refGenomeSeq)
            if (mismatches > maxMismatches)
            {
                sLogger.error(
                    "validation failed: seq({}) and ref genome seq({} of {}) sequence mismatch({}) > $maxMismatches",
                    seq,
                    refGenomeSeq,
                    genomicLocation,
                    mismatches
                )
                return false
            }
            return true
        }

        private fun getChromosomeForRefGenome(chromosome: String, refGenome: RefGenomeSource): String
        {
            return if (refGenome.refGenomeFile().index.hasIndexEntry(chromosome))
                chromosome
            else
                // maybe need to try removing chr
                chromosome.replace("chr", "")
        }

        private fun countMismatches(seq1: String, seq2: String): Int
        {
            require(seq1.length == seq2.length) { "Sequences must be same length" }
            return seq1.zip(seq2).count { it.first != it.second }
        }
    }
}

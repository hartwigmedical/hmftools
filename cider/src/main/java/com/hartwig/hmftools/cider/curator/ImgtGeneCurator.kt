package com.hartwig.hmftools.cider.curator

import com.beust.jcommander.JCommander
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.UnixStyleUsageFormatter
import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.ALIGNMENT_MAX_MISMATCH
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IGKDEL_SEQ
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IGKINTR_SEQ
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IMGT_ANCHOR_LENGTH
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IMGT_V_ANCHOR_INDEX
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.SPECIES
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.jAnchorSignatures
import com.hartwig.hmftools.cider.genes.Contig
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.cider.genes.anchorGenomicLocation
import com.hartwig.hmftools.common.bwa.BwaUtils
import com.hartwig.hmftools.common.cider.IgTcrFunctionality
import com.hartwig.hmftools.common.cider.IgTcrGene
import com.hartwig.hmftools.common.cider.IgTcrGeneFile
import com.hartwig.hmftools.common.cider.IgTcrRegion
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache
import com.hartwig.hmftools.common.gene.GeneData
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.region.BaseRegion
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import htsjdk.samtools.reference.FastaSequenceFile
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import org.apache.logging.log4j.core.config.Configurator
import java.io.File
import kotlin.system.exitProcess

const val ANCHOR_DNA_LENGTH: Int = 30

// This utility processes IMGT IG/TCR genes into format CIDER can use.
// It uses BWA-MEM combined with ensembl to find the genomic location of each IG/TCR genes.
// Download the IMGT sequences from the IMGT website.
// Format for the fasta line is F+ORF+in-frame P nucleotide sequences with IMGT gaps
// See IMGT doc: https://www.imgt.org/genedb/doc
class ImgtGeneCurator
{
    @Parameter(names = ["-imgt"], required = true, description = "Input IMGT fasta file")
    lateinit var inputImgtFasta: String

    @Parameter(names = ["-" + EnsemblDataCache.ENSEMBL_DATA_DIR], required = true, description = EnsemblDataCache.ENSEMBL_DATA_DIR_CFG)
    lateinit var ensemblDataDir: String

    @Parameter(names = ["-ref_genome"], required = true, description = "Reference genome fasta file")
    lateinit var refGenomePath: String

    @Parameter(names = ["-bwa_index_image"], required = true, description = "BWA-MEM index image for the reference genome")
    lateinit var bwaIndexImagePath: String

    @Parameter(names = ["-${BwaUtils.BWA_LIB_PATH}"], required = false)
    var bwaLibPath: String? = null

    @Parameter(names = ["-output"], required = true, description = "Output TSV file")
    lateinit var output: String

    @Parameter(names = ["-threads"], description = "Number of threads")
    var threadCount = 1

    @Parameter(names = ["-log_level"], required = false)
    var logLevel: String = "info"

    data class ImgtGene(val geneName: String, val allele: String, val species: String, val functionality: IgTcrFunctionality,
                        val region: IgTcrRegion?, val sequenceWithGaps: String, val partial: Boolean = false)
    {
        val sequenceWithoutGaps: String get() { return sequenceWithGaps.replace(".", "") }
    }

    fun run(): Int
    {
        val refGenome = RefGenomeSource.loadRefGenome(refGenomePath)
        val refGenomeVersion = deriveRefGenomeVersion(refGenome)

        val ensemblDataCache = EnsemblDataCache(ensemblDataDir, refGenomeVersion)
        val ensemblLoadOk = ensemblDataCache.load(true)

        if (!ensemblLoadOk)
        {
            sLogger.error("Ensembl data cache load failed")
            throw RuntimeException("Ensembl data cache load failed")
        }

        BwaUtils.loadAlignerLibrary(bwaLibPath)

        val imgtGeneList = readImgtGenesFromFasta(inputImgtFasta)
        imgtGeneList.addAll(getCustomGenes())

        val refGenomeDict = "$refGenomePath.dict"
        val geneLocations = findGeneLocations(imgtGeneList, refGenomeDict, bwaIndexImagePath, refGenomeVersion, threadCount, ensemblDataCache)
        require(geneLocations.size == imgtGeneList.size)

        val igTcrGeneList = ArrayList<IgTcrGene>()
        for ((geneData, geneLocation) in imgtGeneList.zip(geneLocations))
        {
            processImgtGene(geneData, geneLocation)?.let { igTcrGeneList.add(it) }
        }

        // validate anchor locations
        sLogger.info("validating anchor locations")
        validateAnchorLocations(igTcrGeneList, refGenomePath)

        IgTcrGeneFile.write(output, igTcrGeneList)

        return 0
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(ImgtGeneCurator::class.java)

        @JvmStatic
        fun main(args: Array<String>)
        {
            // here we have some voodoo to work out if we are being used in pipeline mode or the standalone mode
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
                Configurator.setRootLevel(Level.valueOf(app.logLevel))
                exitProcess(app.run())
            } catch (paramException: ParameterException)
            {
                println("${paramException.message}")
                commander.usage()
                exitProcess(1)
            }
        }

        // https://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html
        fun functionalityFromImgtCode(code: String): IgTcrFunctionality
        {
            return when (code)
            {
                "F" -> IgTcrFunctionality.FUNCTIONAL
                "ORF" -> IgTcrFunctionality.ORF
                "P" -> IgTcrFunctionality.PSEUDOGENE
                else -> throw IllegalArgumentException("unrecognised Imgt functionality: $code")
            }
        }

        // IGH has three (alpha, delta and gamma) or four (epsilon and mu) constant domains (CH1 to CH4)
        fun igTcrRegionFromImgtCode(region: String): IgTcrRegion?
        {
            return when (region)
            {
                "V-REGION" -> IgTcrRegion.V_REGION
                "D-REGION" -> IgTcrRegion.D_REGION
                "J-REGION" -> IgTcrRegion.J_REGION
                "C-REGION" -> IgTcrRegion.CONSTANT // IGK / IGL
                "CH1" -> IgTcrRegion.CONSTANT // IGH
                "EX1" -> IgTcrRegion.CONSTANT // TCR
                else -> null
            }
        }

        /*
        The FASTA header of IMGT/GENE-DB reference sequences is standardized. It contains 15 fields separated by '|':

        1. IMGT/LIGM-DB accession number(s)
        2. IMGT gene and allele name
        3. species
        4. IMGT allele functionality
        5. exon(s), region name(s), or extracted label(s)
        6. start and end positions in the IMGT/LIGM-DB accession number(s)
        7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
        8. codon start, or 'NR' (not relevant) for non coding labels
        9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
        10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
        11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
        12. number of amino acids (AA): this field indicates that the sequence is in amino acids
        13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
        14. partial (if it is)
        15. reverse complementary (if it is)
         */
        // >IMGT000128|IGHA1*06|Homo sapiens|F|M|g,1187575..1187786|213 nt|1|+1| | | |213+0=213| | |
        // https://www.imgt.org/IMGTScientificChart/SequenceDescription/IMGTfunctionality.html
        fun parseImgtGeneData(seqName: String, sequenceWithGaps: String): ImgtGene
        {
            val tokens = seqName.split('|')
            require(tokens.size > 10)

            val geneAllele = tokens[1].split('*')
            require(geneAllele.size == 2)

            // we don't distinguish between F, (F) and [F]
            val functionality = functionalityFromImgtCode(tokens[3].replace(Regex("[()\\[\\]]"), ""))

            val partial = tokens[13].contains("partial")

            return ImgtGene(
                geneName = geneAllele[0], allele = geneAllele[1], species = tokens[2], functionality = functionality,
                region = igTcrRegionFromImgtCode(tokens[4]), sequenceWithGaps = sequenceWithGaps, partial = partial)
        }

        fun readImgtGenesFromFasta(imgtFastaPath: String): ArrayList<ImgtGene>
        {
            val imgtFastaFile = FastaSequenceFile(File(imgtFastaPath), false)

            val imgtGeneList = ArrayList<ImgtGene>()

            while (true)
            {
                val sequence = imgtFastaFile.nextSequence()

                if (sequence == null)
                    break

                val imgtGeneData = parseImgtGeneData(sequence.name, sequence.baseString.uppercase())
                imgtGeneList.add(imgtGeneData)

                sLogger.info("imgt gene: {}", imgtGeneData)
            }

            return imgtGeneList
        }

        fun getCustomGenes(): List<ImgtGene> =
            // add IGKINTR and IGKDEL
            listOf(
                ImgtGene(
                    geneName = VJGeneType.IGKINTR, allele = "01", species = SPECIES, functionality = IgTcrFunctionality.ORF,
                    region = IgTcrRegion.V_REGION, sequenceWithGaps = IGKINTR_SEQ
                ),
                ImgtGene(
                    geneName = VJGeneType.IGKDEL, allele = "01", species = SPECIES, functionality = IgTcrFunctionality.ORF,
                    region = IgTcrRegion.J_REGION, sequenceWithGaps = IGKDEL_SEQ
                ))

        fun processImgtGene(geneData: ImgtGene, geneLocation: GenomicLocation?): IgTcrGene?
        {
            // we filter by species
            if (geneData.species != SPECIES)
            {
                return null
            }

            // get the gene region
            val region = geneData.region ?: return null

            val anchorData: Pair<String, GenomicLocation?>? = when (region)
            {
                IgTcrRegion.V_REGION -> findVAnchor(geneData, geneLocation)
                IgTcrRegion.J_REGION -> findJAnchor(geneData, geneLocation)
                else -> null
            }

            require(anchorData == null || geneLocation?.contig == anchorData.second?.contig)
            return IgTcrGene(
                geneData.geneName,
                geneData.allele,
                region,
                geneData.functionality,
                geneLocation?.contig?.name,
                geneLocation?.position,
                geneLocation?.strand,
                anchorData?.first,
                anchorData?.second?.position
            )
        }

        fun findVAnchor(geneData: ImgtGene, geneLocation: GenomicLocation?): Pair<String, GenomicLocation?>?
        {
            // some sequences are longer
            val seqWithGaps = geneData.sequenceWithGaps

            if (seqWithGaps.length < IMGT_V_ANCHOR_INDEX)
            {
                sLogger.log(
                    if (geneData.functionality == IgTcrFunctionality.FUNCTIONAL && !geneData.partial) Level.ERROR else Level.INFO,
                    "Cannot find V anchor, sequence too short. {}", geneData)
                return null
            }

            var anchor = seqWithGaps.substring(IMGT_V_ANCHOR_INDEX, Math.min(IMGT_V_ANCHOR_INDEX + IMGT_ANCHOR_LENGTH, seqWithGaps.length))

            // if anchor is too short we remove
            if (anchor.length < IMGT_ANCHOR_LENGTH)
            {
                // skip this one
                sLogger.log(
                    if (geneData.functionality == IgTcrFunctionality.FUNCTIONAL && !geneData.partial) Level.ERROR else Level.INFO,
                    "V anchor:{} too short for {}", anchor, geneData
                )
                return null
            }

            // anchor cannot contain .
            if (anchor.contains("."))
            {
                sLogger.error("V anchor: {} contains \".\", gene: {}", anchor, geneData)
                return null
            }

            // now we find the anchor index again using the one without gap
            val anchorIndex = geneData.sequenceWithoutGaps.indexOf(anchor)

            // if the anchor is < 30 bases long, we fill it in with the last C which is TGT
            // it is missing in some
            // this is hacky
            if (anchor.length < ANCHOR_DNA_LENGTH)
            {
                sLogger.warn("V anchor for {} too short, adding TGT to end", geneData)
                anchor = anchor.take(27) + "TGT"
            }

            // for v gene the anchor is near the end, due to introns we much find the sequence from the end backwards
            val anchorOffsetFromEnd = geneData.sequenceWithoutGaps.length - anchorIndex - anchor.length

            var anchorLocation: GenomicLocation? = null

            if (geneLocation != null && anchorIndex != -1)
            {
                anchorLocation = if (geneLocation.strand == Strand.FORWARD)
                {
                    geneLocation.copy(
                        position = BaseRegion(
                            geneLocation.position.end() - anchorOffsetFromEnd - anchor.length + 1,
                            geneLocation.position.end() - anchorOffsetFromEnd
                        ))
                } else
                {
                    geneLocation.copy(
                        position = BaseRegion(
                            geneLocation.position.start() + anchorOffsetFromEnd,
                            geneLocation.position.start() + anchorOffsetFromEnd + anchor.length - 1
                        ))
                }
            }

            val aaSeq = Codons.aminoAcidFromBases(anchor)

            // v gene
            sLogger.info(
                "V gene: {}, anchor: {}, offset from end: {}, anchor AA: {}",
                geneData.geneName,
                anchor,
                anchorOffsetFromEnd,
                aaSeq
            )

            return Pair(anchor, anchorLocation)
        }

        fun findJAnchor(geneData: ImgtGene, geneLocation: GenomicLocation?): Pair<String, GenomicLocation?>?
        {
            // J gene rules
            // 30 base sequence starting with TGGGG (W) or TTTG and TTCG (F)
            val seqWithGaps = geneData.sequenceWithGaps
            val anchorIndex: Int = seqWithGaps.indexOfAny(jAnchorSignatures(geneData.geneName))

            if (anchorIndex <= 0)
            {
                sLogger.log(
                    if (geneData.functionality == IgTcrFunctionality.FUNCTIONAL && !geneData.partial) Level.ERROR else Level.INFO,
                    "J gene: {} cannot find anchor", geneData)
                return null
            }

            val anchor = seqWithGaps.substring(anchorIndex, Math.min(anchorIndex + IMGT_ANCHOR_LENGTH, seqWithGaps.length))

            // cannot contain .
            if (anchor.contains("."))
            {
                sLogger.error("J gene: {}, anchor({}) contains .", geneData, anchor)
                return null
            }

            var anchorLocation: GenomicLocation? = null

            if (geneLocation != null)
            {
                anchorLocation = if (geneLocation.strand == Strand.FORWARD)
                {
                    geneLocation.copy(
                        position = BaseRegion(
                            geneLocation.position.start() + anchorIndex,
                            geneLocation.position.start() + anchorIndex + anchor.length - 1
                        ))
                } else
                {
                    geneLocation.copy(
                        position = BaseRegion(
                            geneLocation.position.end() - anchorIndex - anchor.length + 1,
                            geneLocation.position.end() - anchorIndex
                        ))
                }
            }

            val aaSeq = Codons.aminoAcidFromBases(anchor)

            sLogger.info("J gene: {}, anchor: {}, offset from start: {}, anchor AA: {}", geneData.geneName, anchor, anchorIndex, aaSeq)

            return Pair(anchor, anchorLocation)
        }

        fun findGeneLocations(genes: List<ImgtGene>, refGenomeDict: String, bwaIndexImage: String,
                              refGenomeVersion: RefGenomeVersion, numThreads: Int, ensemblDataCache: EnsemblDataCache)
            : List<GenomicLocation?>
        {
            // First, apply any hardcoded location overrides, which take precedence.
            val locations = genes.map { ImgtGeneCuratorSettings.getGenomicLocationOverride(it.geneName, refGenomeVersion) }
                .toMutableList()

            // Find locations of constant genes.
            for ((i, gene) in genes.withIndex())
            {
                if (locations[i] == null && gene.region == IgTcrRegion.CONSTANT)
                {
                    val ensemblGene = ensemblDataCache.getGeneDataByName(gene.geneName)
                    if (ensemblGene != null)
                    {
                        locations[i] = toGenomicLocation(ensemblGene, refGenomeVersion)
                    }
                }
            }

            val nonConstantGenes = genes.withIndex().filter { locations[it.index] == null && it.value.region != IgTcrRegion.CONSTANT }
            // V / D / J gene use alignment, constant genes use ensembl
            // reason is that ensembl is easier, and we do not need the alt locations for the constant genes. Another reason
            // is that we do not need to be very precise with the location of the anchor for constant genes
            // if we use ensembl for V / J gene, we need to use the fasta file to validate the anchor location is precise
            val nonConstantGeneLocations = alignForGenomicLocation(
                nonConstantGenes.map { it.value }, refGenomeDict, bwaIndexImage, refGenomeVersion, numThreads, ensemblDataCache)
            for ((gene, location) in nonConstantGenes.zip(nonConstantGeneLocations))
            {
                locations[gene.index] = location
            }

            return locations
        }

        // for each imgt gene segment, use alignment to find the genomic location
        fun alignForGenomicLocation(imgtGeneList: List<ImgtGene>, refGenomeDict: String, bwaIndexImage: String,
                                    refGenomeVersion: RefGenomeVersion, numThreads: Int, ensemblDataCache: EnsemblDataCache)
            : List<GenomicLocation?>
        {
            val locations = ArrayList<GenomicLocation?>()

            val alignments = AlignmentUtil.runBwaMem(
                imgtGeneList.map { it.sequenceWithoutGaps },
                refGenomeDict, bwaIndexImage,
                0, numThreads)

            for ((geneData, geneAlignments) in imgtGeneList.zip(alignments))
            {
                sLogger.debug("Processing alignments for gene {}", geneData)

                val alignmentsWithLocations = geneAlignments
                    .filter { it.editDistance <= ALIGNMENT_MAX_MISMATCH &&
                            it.queryAlignEnd - it.queryAlignStart + 1 >= it.querySeq.length - ALIGNMENT_MAX_MISMATCH }
                    .mapNotNull {AlignmentUtil.toGenomicLocation(it)?.let { location -> Pair(it, location) } }
                    // Order by alignment score then primary assembly first, with position as a tie-breaker
                    .sortedWith(compareBy(
                        { -it.first.alignmentScore }, { !it.second.inPrimaryAssembly }, { it.first.refContig }, { it.first.refStart }))

                val ensemblGene = ensemblDataCache.getGeneDataByName(toEnsemblGeneName(geneData.geneName))
                var bestAlignment: AlignmentUtil.BwaMemAlignment? = null

                // TODO: should resolve equal alignments similar to how it's done in cider main app

                for ((alignment, location) in alignmentsWithLocations)
                {
                    sLogger.debug("considering alignment: {}", alignment)

                    if (bestAlignment == null) {
                        bestAlignment = alignment
                    }

                    else if (bestAlignment.alignmentScore == alignment.alignmentScore)
                    {
                        // if they have same score, we want to choose the one that overlaps with ensembl gene
                        sLogger.debug("considering ensembl gene {}", ensemblGene)
                        if (ensemblGene != null &&
                            // TODO: check if it's chromosome or contig
                            ensemblGene.Chromosome == location.contig.name &&
                            ensemblGene.forwardStrand() == (location.strand == Strand.FORWARD) &&
                            ensemblGene.GeneStart <= location.position.end() &&
                            ensemblGene.GeneEnd >= location.position.start())
                        {
                            bestAlignment = alignment
                            sLogger.debug("gene: {}*{}, using alignment that overlaps with ensembl", geneData.geneName, geneData.allele)
                        }
                    }
                }

                var location: GenomicLocation? = null
                if (bestAlignment == null)
                {
                    if (ensemblGene != null)
                    {
                        sLogger.error("gene: {}*{}, no full match yet has ensembl", geneData.geneName, geneData.allele)

                        if (geneData.region == IgTcrRegion.D_REGION)
                        {
                            // use ensembl for D region
                            location = toGenomicLocation(ensemblGene, refGenomeVersion)
                        }
                    }
                    sLogger.info("gene: {}*{}, no full match", geneData.geneName, geneData.allele)
                }
                else
                {
                    location = alignmentToQueryGenomicLocation(bestAlignment)
                    sLogger.info("gene: {}*{}, match: {}, gene loc: {}", geneData.geneName, geneData.allele, bestAlignment, location)
                }
                locations.add(location)
            }

            return locations
        }

        fun toEnsemblGeneName(name: String) = name.replace("/", "").uppercase()

        fun toGenomicLocation(ensemblGene: GeneData, refGenomeVersion: RefGenomeVersion) : GenomicLocation
        {
            // TODO: what is ensembleGene.Chromosome compared to contig?
            val contig = Contig(refGenomeVersion.versionedChromosome(ensemblGene.Chromosome))
            return GenomicLocation(
                contig, BaseRegion(ensemblGene.GeneStart,ensemblGene.GeneEnd),
                Strand.valueOf(ensemblGene.Strand.toInt()))
        }

        fun alignmentToQueryGenomicLocation(alignment: AlignmentUtil.BwaMemAlignment) : GenomicLocation
        {
            // we need to correct for the ends to make sure things align properly
            val startExtend = alignment.queryAlignStart - 1
            val endExtend = alignment.querySeq.length - alignment.queryAlignEnd
            val alignGenomicLoc = AlignmentUtil.toGenomicLocation(alignment)!!

            if (startExtend == 0 && endExtend == 0)
            {
                return alignGenomicLoc
            }

            return if (alignGenomicLoc.strand == Strand.FORWARD)
            {
                alignGenomicLoc.copy(
                    position = BaseRegion(
                        alignGenomicLoc.position.start() - startExtend,
                        alignGenomicLoc.position.end() + endExtend
                    ))
            }
            else
            {
                alignGenomicLoc.copy(
                    position = BaseRegion(
                        alignGenomicLoc.position.start() - endExtend,
                        alignGenomicLoc.position.end() + startExtend
                    ))
            }
        }

        fun validateAnchorLocations(igTcrGeneList: List<IgTcrGene>, refGenomeFastaPath: String)
        {
            val refGenomeFile = IndexedFastaSequenceFile(File(refGenomeFastaPath))
            val genomicLocationValidator = GenomicLocationValidator(refGenomeFile)

            for (gene in igTcrGeneList)
            {
                val anchorGenomicLocation = gene.anchorGenomicLocation()
                if (anchorGenomicLocation != null && anchorGenomicLocation.inPrimaryAssembly)
                {
                    if (!genomicLocationValidator.validateAgainstRefGenome(gene.anchorSequence!!, anchorGenomicLocation))
                    {
                        sLogger.error("gene: {} anchor location: {} does not match anchor seq: {}",
                            gene.geneAllele, gene.anchorPosition, gene.anchorSequence)
                        throw RuntimeException()
                    }
                }
            }
        }
    }
}

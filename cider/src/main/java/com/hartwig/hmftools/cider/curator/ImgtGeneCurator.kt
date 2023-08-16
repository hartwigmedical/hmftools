package com.hartwig.hmftools.cider.curator

import com.beust.jcommander.JCommander
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.UnixStyleUsageFormatter
import com.google.common.collect.Multimap
import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.cider.CiderConstants.BLAST_REF_GENOME_VERSION
import com.hartwig.hmftools.cider.VJGeneType.Companion.IGKDEL
import com.hartwig.hmftools.cider.blastn.BlastnMatch
import com.hartwig.hmftools.cider.blastn.BlastnMatch.Companion.PRIMARY_ASSEMBLY_NAME
import com.hartwig.hmftools.cider.blastn.BlastnRunner
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.cider.genes.IgTcrGeneFile
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.utils.Doubles
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import htsjdk.samtools.reference.FastaSequenceFile
import org.apache.logging.log4j.Level
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.Comparator
import kotlin.system.exitProcess

const val ANCHOR_DNA_LENGTH: Int = 30

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
// -outputV38 igtcr_gene.38.tsv
// -outputV37 igtcr_gene.37.tsv
//
// Format for the fasta line is F+ORF+in-frame P nucleotide sequences with IMGT gaps
// See IMGT doc: https://www.imgt.org/genedb/doc
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

    @Parameter(names = ["-output_v38"], required = true, description = "Output TSV file for HG38")
    lateinit var outputV38: String

    @Parameter(names = ["-output_v37"], description = "Output TSV file for HG37")
    var outputV37: String? = null

    @Parameter(names = ["-threads"], description = "Number of threads")
    var threadCount = CiderParams.DEFAULT_THREADS

    @Parameter(names = ["-workdir"], description = "Number of threads")
    lateinit var workdir: String

    data class ImgtGeneData(val geneName: String, val allele: String, val species: String, val functionality: IgTcrFunctionality,
                            val region: String, val sequenceWithGaps: String, val partial: Boolean = false,
                            var genomicLocation: GenomicLocation? = null, var blastnMatch: BlastnMatch? = null)
    {
        val sequenceWithoutGaps: String get() { return sequenceWithGaps.replace(".", "") }
    }

    fun run(): Int
    {
        val ensemblDataCache = EnsemblDataCache(ensemblDataDir, BLAST_REF_GENOME_VERSION)
        val ensemblLoadOk = ensemblDataCache.load(true)

        if (!ensemblLoadOk)
        {
            sLogger.error("Ensembl data cache load failed")
            throw RuntimeException("Ensembl data cache load failed")
        }

        val imgtGeneDataList: List<ImgtGeneData> = readGeneDataFromFasta(inputImgtFasta)

        blastForGenomicLocation(imgtGeneDataList, blast, blastDb, threadCount, workdir, ensemblDataCache)

        val igTcrGeneList = ArrayList<IgTcrGene>()

        for (geneData in imgtGeneDataList)
        {
            processImgtGeneData(geneData)?.let { igTcrGeneList.add(it) }
        }

        IgTcrGeneFile.write(outputV38, igTcrGeneList)

        if (outputV37 != null)
        {
            val igTcrGeneListV37 = convertToV37(igTcrGeneList)
            IgTcrGeneFile.write(outputV37!!, igTcrGeneListV37)
        }

        return 0
    }

    companion object
    {
        const val SPECIES = "Homo sapiens"
        const val IMGT_V_ANCHOR_INDEX = 282
        const val IMGT_ANCHOR_LENGTH = 30
        const val BLASTN_EVALUE_CUTOFF = 1000.0
        const val BLASTN_MAX_MISMATCH = 1

        // shift the V anchor index such that it starts at the first base
        val IGKINTR_SEQ = ".".repeat(IMGT_V_ANCHOR_INDEX) + "CACCGCGCTCTTGGGGCAGCCGCCTTGCCGCTAGTGGCCGTGGCCACCCTGTGTCTGCCCGATT"
        val IGKDEL_SEQ = "GGAGCCCTAGTGGCAGCCCAGGGCGACTCCTCATGAGTCTGCAGCTGCATTTTTGCCATATCCACTATTTGGAGTCTGACCTCCCTAGGAAGCCTCCCTGC"

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
                exitProcess(app.run())
            } catch (paramException: ParameterException)
            {
                println("${paramException.message}")
                commander.usage()
                exitProcess(1)
            }
        }

        fun jAnchorSignatures(geneName: String) : List<String>
        {
            return if (geneName.startsWith("IGHJ"))
                listOf("TGGGG")
            else if (geneName == IGKDEL)
                listOf("GCCC")
            else
                listOf("TTTG", "TTCG")
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
        fun parseGeneData(seqName: String, sequenceWithGaps: String): ImgtGeneData
        {
            val tokens = seqName.split('|')
            require(tokens.size > 10)

            val geneAllele = tokens[1].split('*')
            require(geneAllele.size == 2)

            // we don't distinguish between F, (F) and [F]
            val functionality = functionalityFromImgtCode(tokens[3].replace(Regex("[()\\[\\]]"), ""))

            val partial = tokens[13].contains("partial")

            return ImgtGeneData(
                geneName = geneAllele[0], allele = geneAllele[1], species = tokens[2], functionality = functionality,
                region = tokens[4], sequenceWithGaps = sequenceWithGaps, partial = partial)
        }

        fun readGeneDataFromFasta(imgtFastaPath: String): List<ImgtGeneData>
        {
            val imgtFastaFile = FastaSequenceFile(File(imgtFastaPath), false)

            val imgtGeneDataList: MutableList<ImgtGeneData> = ArrayList()

            while (true)
            {
                val sequence = imgtFastaFile.nextSequence()

                if (sequence == null)
                    break

                val imgtGeneData = parseGeneData(sequence.name, sequence.baseString.uppercase())
                imgtGeneDataList.add(imgtGeneData)

                sLogger.info("imgt gene: {}", imgtGeneData)
            }

            // add IGKINTR and IGKDEL
            imgtGeneDataList.add(ImgtGeneData(
                geneName = VJGeneType.IGKINTR, allele = "01", species = SPECIES, functionality = IgTcrFunctionality.FUNCTIONAL,
                region = "V-REGION", sequenceWithGaps = IGKINTR_SEQ
            ))

            imgtGeneDataList.add(ImgtGeneData(
                geneName = VJGeneType.IGKDEL, allele = "01", species = SPECIES, functionality = IgTcrFunctionality.FUNCTIONAL,
                region = "J-REGION", sequenceWithGaps = IGKDEL_SEQ
            ))

            return imgtGeneDataList
        }

        fun processImgtGeneData(geneData: ImgtGeneData): IgTcrGene?
        {
            // we filter by species
            if (geneData.species != SPECIES)
            {
                return null
            }

            // get the gene region
            val region = igTcrRegionFromImgtCode(geneData.region) ?: return null

            val anchorData: Pair<String, GenomicLocation?>? = when (region)
            {
                IgTcrRegion.V_REGION -> findVAnchor(geneData)
                IgTcrRegion.J_REGION -> findJAnchor(geneData)
                else -> null
            }

            return IgTcrGene(
                geneData.geneName,
                geneData.allele,
                region,
                geneData.functionality,
                geneData.genomicLocation,
                anchorData?.first,
                anchorData?.second
            )
        }

        fun findVAnchor(geneData: ImgtGeneData): Pair<String, GenomicLocation?>?
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

            val geneLocation = geneData.genomicLocation

            if (geneLocation != null && anchorIndex != -1)
            {
                anchorLocation = if (geneLocation.strand == Strand.FORWARD)
                {
                    geneLocation.copy(
                        posStart = geneLocation.posEnd - anchorOffsetFromEnd - anchor.length + 1,
                        posEnd = geneLocation.posEnd - anchorOffsetFromEnd
                    )
                } else
                {
                    geneLocation.copy(
                        posStart = geneLocation.posStart + anchorOffsetFromEnd,
                        posEnd = geneLocation.posStart + anchorOffsetFromEnd + anchor.length - 1
                    )
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

        fun findJAnchor(geneData: ImgtGeneData): Pair<String, GenomicLocation?>?
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

            val geneLocation = geneData.genomicLocation

            if (geneLocation != null)
            {
                anchorLocation = if (geneLocation.strand == Strand.FORWARD)
                {
                    geneLocation.copy(
                        posStart = geneLocation.posStart + anchorIndex,
                        posEnd = geneLocation.posStart + anchorIndex + anchor.length - 1
                    )
                } else
                {
                    geneLocation.copy(
                        posStart = geneLocation.posEnd - anchorIndex - anchor.length + 1,
                        posEnd = geneLocation.posEnd - anchorIndex
                    )
                }
            }

            val aaSeq = Codons.aminoAcidFromBases(anchor)

            sLogger.info("J gene: {}, anchor: {}, offset from start: {}, anchor AA: {}", geneData.geneName, anchor, anchorIndex, aaSeq)

            return Pair(anchor, anchorLocation)
        }

        // for each imgt gene segment, use blastn to find the genomic location
        fun blastForGenomicLocation(imgtGeneDataList: List<ImgtGeneData>, blastn: String, blastDb: String, numThreads: Int, workdir: String,
                                    ensemblDataCache: EnsemblDataCache)
        {
            // assign a key to each VDJ, such that we can keep track of them
            var key = 0
            val keyToGeneDataMap: Map<Int, ImgtGeneData> = imgtGeneDataList.associateBy { ++key }

            val blastnResults: Multimap<Int, BlastnMatch> = BlastnRunner.runBlastn(
                "igmt", blastn, blastDb,
                keyToGeneDataMap.mapValues { geneData -> geneData.value.sequenceWithoutGaps },
                workdir, numThreads, BLASTN_EVALUE_CUTOFF,
                true
            )

            // process the blastnResults
            for ((k, geneData) in keyToGeneDataMap)
            {
                // Filter by full match and not too many mismatches, then
                // We sort the matches by
                // 1. match quality
                // 2. primary assembly
                val matches = blastnResults[k]
                    .filter { m: BlastnMatch -> m.numMismatch <= BLASTN_MAX_MISMATCH && m.alignmentLength == m.querySeqLen }
                    .sortedWith(Comparator.comparingDouble { m: BlastnMatch -> m.expectedValue }
                        .thenComparingInt { m: BlastnMatch -> if (m.subjectTitle.endsWith(PRIMARY_ASSEMBLY_NAME)) 0 else 1 })

                // find the gene in the ensembl
                val ensemblGene = ensemblDataCache.getGeneDataByName(geneData.geneName)

                for (match in matches)
                {
                    val matchLocation = match.toGenomicLocation()

                    if (matchLocation == null)
                    {
                        continue
                    }

                    if (geneData.blastnMatch == null)
                    {
                        geneData.blastnMatch = match
                        geneData.genomicLocation = matchLocation
                    }

                    else if (Doubles.equal(geneData.blastnMatch!!.bitScore, match.bitScore))
                    {
                        // if they have same score, we want to choose the one that overlaps with ensembl gene
                        if (ensemblGene != null &&
                            ensemblGene.Chromosome == matchLocation.chromosome &&
                            ensemblGene.forwardStrand() == (matchLocation.strand == Strand.FORWARD) &&
                            ensemblGene.GeneStart <= matchLocation.posEnd &&
                            ensemblGene.GeneEnd >= matchLocation.posStart)
                        {
                            geneData.blastnMatch = match
                            geneData.genomicLocation = matchLocation
                            sLogger.debug("gene: {}*{}, using match that overlaps with ensembl", geneData.geneName, geneData.allele)
                        }
                    }
                }

                if (geneData.blastnMatch == null)
                {
                    sLogger.info("gene: {}*{}, no full match", geneData.geneName, geneData.allele)
                }
                else
                {
                    sLogger.info("gene: {}*{}, match: {}", geneData.geneName, geneData.allele, geneData.blastnMatch)
                }
            }
        }

        fun convertToV37(igTcrGeneList: List<IgTcrGene>): List<IgTcrGene>
        {
            val genomicLiftOver = GenomeLiftoverCache(true, false)
            val igTcrGeneV37List = ArrayList<IgTcrGene>()

            for (gene in igTcrGeneList)
            {
                val geneLocationV37: GenomicLocation? = if (gene.geneLocation != null)
                    convertGenomicLocationTo37(gene.geneLocation, genomicLiftOver)
                else null

                val anchorLocation: GenomicLocation? = if (gene.anchorLocation != null)
                    convertGenomicLocationTo37(gene.anchorLocation, genomicLiftOver)
                else null

                igTcrGeneV37List.add(gene.copy(geneLocation = geneLocationV37, anchorLocation = anchorLocation))
            }
            return igTcrGeneV37List
        }

        fun convertGenomicLocationTo37(genomicLocation: GenomicLocation, genomicLiftOver: GenomeLiftoverCache): GenomicLocation?
        {
            if (genomicLocation.isPrimaryAssembly)
            {
                val v37PosStart: Int = genomicLiftOver.convertPositionTo37(genomicLocation.chromosome, genomicLocation.posStart)
                val v37PosEnd: Int = genomicLiftOver.convertPositionTo37(genomicLocation.chromosome, genomicLocation.posEnd)

                if (v37PosStart != -1 && v37PosEnd != -1)
                {
                    return genomicLocation.copy(chromosome = RefGenomeVersion.V37.versionedChromosome(genomicLocation.chromosome),
                        posStart = v37PosStart, posEnd = v37PosEnd)
                }
            }
            return null
        }
    }
}

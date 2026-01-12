package com.hartwig.hmftools.cider.curator

import com.beust.jcommander.JCommander
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.UnixStyleUsageFormatter
import com.google.common.collect.Multimap
import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.cider.CiderConstants.BLAST_REF_GENOME_VERSION
import com.hartwig.hmftools.cider.IgTcrGene.Companion.toCommonIgTcrGene
import com.hartwig.hmftools.cider.AlignmentUtil
import com.hartwig.hmftools.cider.AlignmentUtil.BLASTN_PRIMARY_ASSEMBLY_NAME
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.BLASTN_EVALUE_CUTOFF
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.BLASTN_MAX_MISMATCH
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IGKDEL_SEQ
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IGKINTR_SEQ
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IMGT_ANCHOR_LENGTH
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.IMGT_V_ANCHOR_INDEX
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.SPECIES
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.jAnchorSignatures
import com.hartwig.hmftools.cider.curator.ImgtGeneCuratorSettings.liftOverBlacklist
import com.hartwig.hmftools.cider.genes.GenomicLocation
import com.hartwig.hmftools.common.blastn.BlastnMatch
import com.hartwig.hmftools.common.cider.IgTcrGeneFile
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache
import com.hartwig.hmftools.common.gene.GeneData
import com.hartwig.hmftools.common.genome.refgenome.GenomeLiftoverCache
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions.stripChrPrefix
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.utils.Doubles
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import htsjdk.samtools.liftover.LiftOver
import htsjdk.samtools.reference.FastaSequenceFile
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.samtools.util.Interval
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
// -liftover_chain /data/hg38ToHg19.over.chain.gz
// -ref_genome_v38 /data/resources/bucket/reference_genome/38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
// -ref_genome_v37 /data/resources/bucket/reference_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta
// -outputV38 igtcr_gene.38.tsv
// -outputV37 igtcr_gene.37.tsv
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

    @Parameter(names = ["-output_v38"], required = true, description = "Output TSV file for HG38")
    lateinit var outputV38: String

    @Parameter(names = ["-output_v37"], required = true, description = "Output TSV file for HG37")
    lateinit var outputV37: String

    @Parameter(names = ["-threads"], description = "Number of threads")
    var threadCount = 1

    @Parameter(names = ["-workdir"], description = "Number of threads")
    lateinit var workdir: String

    data class ImgtGeneData(val geneName: String, val allele: String, val species: String, val functionality: IgTcrFunctionality,
                            val region: IgTcrRegion?, val sequenceWithGaps: String, val partial: Boolean = false,
                            var genomicLocationV38: GenomicLocation? = null)
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

        findGenomicLocatons(imgtGeneDataList, blast, blastDb, threadCount, workdir, ensemblDataCache)

        //val genomicLiftOver = GenomeLiftoverCache(true, false)
        val genomicLiftOver = LiftOver(File(liftOverChainFile))

        val igTcrGeneListV38 = ArrayList<IgTcrGene>()

        for (geneData in imgtGeneDataList)
        {
            processImgtGeneData(geneData)?.let { igTcrGeneListV38.add(it) }
        }

        // now convert to V37
        val igTcrGeneListV37 = igTcrGeneListV38.map { gene -> convertGeneTo37(gene, genomicLiftOver) }

        // validate anchor locations
        sLogger.info("validating V38 anchor locations")
        validateAnchorLocations(igTcrGeneListV38, refGenomeV38)
        sLogger.info("validating V37 anchor locations")
        validateAnchorLocations(igTcrGeneListV37, refGenomeV37)

        IgTcrGeneFile.write(outputV38, igTcrGeneListV38.map { o -> toCommonIgTcrGene(o) })
        IgTcrGeneFile.write(outputV37, igTcrGeneListV37.map { o -> toCommonIgTcrGene(o) })

        return 0
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(ImgtGeneCurator::class.java)
        val genomicLiftOverHmf = GenomeLiftoverCache(true)

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
                region = igTcrRegionFromImgtCode(tokens[4]), sequenceWithGaps = sequenceWithGaps, partial = partial)
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
                if (imgtGeneData.species != SPECIES)
                {
                    continue
                }
                if (imgtGeneData.region == null) {
                    continue
                }
                imgtGeneDataList.add(imgtGeneData)

                sLogger.info("imgt gene: {}", imgtGeneData)
            }

            // add IGKINTR and IGKDEL
            imgtGeneDataList.add(ImgtGeneData(
                geneName = VJGeneType.IGKINTR, allele = "01", species = SPECIES, functionality = IgTcrFunctionality.ORF,
                region = IgTcrRegion.V_REGION, sequenceWithGaps = IGKINTR_SEQ
            ))

            imgtGeneDataList.add(ImgtGeneData(
                geneName = VJGeneType.IGKDEL, allele = "01", species = SPECIES, functionality = IgTcrFunctionality.ORF,
                region = IgTcrRegion.J_REGION, sequenceWithGaps = IGKDEL_SEQ
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
            val region = geneData.region ?: return null

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
                geneData.genomicLocationV38,
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

            val geneLocation = geneData.genomicLocationV38

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

            val geneLocation = geneData.genomicLocationV38

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

        fun findGenomicLocatons(imgtGeneDataList: List<ImgtGeneData>, blastn: String, blastDb: String, numThreads: Int, workdir: String,
                                 ensemblDataCache: EnsemblDataCache)
        {
            val constantGenes = imgtGeneDataList.filter { gene -> gene.region == IgTcrRegion.CONSTANT }
            val nonConstantGenes = imgtGeneDataList.filter { gene -> gene.region != IgTcrRegion.CONSTANT }

            // V / D / J gene use blast, constant genes use ensembl
            // reason is that ensembl is easier, and we do not need the alt locations for the constant genes. Another reason
            // is that we do not need to be very precise with the location of the anchor for constant genes
            // if we use ensembl for V / J gene, we need to use the fasta file to validate the anchor location is precise
            blastForGenomicLocation(nonConstantGenes, blastn, blastDb, numThreads, workdir, ensemblDataCache)

            for (geneData in constantGenes)
            {
                val ensemblGene = ensemblDataCache.getGeneDataByName(geneData.geneName) ?: continue
                geneData.genomicLocationV38 = toGenomicLocation(ensemblGene)
            }

            // also apply genomic location overrides
            for (geneData in imgtGeneDataList)
            {
                if (geneData.genomicLocationV38 == null)
                {
                    val geneLocationOverride = ImgtGeneCuratorSettings.getGenomicLocationOverrides(geneData.geneName)
                    if (geneLocationOverride != null)
                    {
                        geneData.genomicLocationV38 = geneLocationOverride
                    }
                }
            }
        }

        // for each imgt gene segment, use blastn to find the genomic location
        fun blastForGenomicLocation(imgtGeneDataList: List<ImgtGeneData>, blastn: String, blastDb: String, numThreads: Int, workdir: String,
                                    ensemblDataCache: EnsemblDataCache)
        {
            // assign a key to each VDJ, such that we can keep track of them
            var key = 0
            val keyToGeneDataMap: Map<Int, ImgtGeneData> = imgtGeneDataList.associateBy { ++key }

            val blastnResults: Multimap<Int, BlastnMatch> = AlignmentUtil.runBlastn(
                "imgt", blastn, blastDb,
                keyToGeneDataMap.mapValues { geneData -> geneData.value.sequenceWithoutGaps },
                workdir, numThreads, BLASTN_EVALUE_CUTOFF)

            println("GeneBlastnData:Name\tAllele\tRegion\tFunctionality\tSequence\tBlastnId\tCiderLocation")

            // process the blastnResults
            for ((k, geneData) in keyToGeneDataMap)
            {
                // Filter by full match and not too many mismatches, then
                // We sort the matches by
                // 1. match quality
                // 2. primary assembly
                val matches = blastnResults[k]
                    .filter { m: BlastnMatch -> m.numMismatch <= BLASTN_MAX_MISMATCH &&
                            m.alignmentLength >= (m.querySeqLen - BLASTN_MAX_MISMATCH) }
                    .sortedWith(Comparator.comparingDouble { m: BlastnMatch -> m.expectedValue }
                        .thenComparingInt { m: BlastnMatch -> if (m.subjectTitle.endsWith(BLASTN_PRIMARY_ASSEMBLY_NAME)) 0 else 1 })

                // find the gene in the ensembl
                val ensemblGene = ensemblDataCache.getGeneDataByName(geneData.geneName)
                var bestMatch: BlastnMatch? = null

                for (match in matches)
                {
                    val matchLocation = AlignmentUtil.toGenomicLocation(match)

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
                            sLogger.debug("gene: {}*{}, using match that overlaps with ensembl", geneData.geneName, geneData.allele)
                        }
                    }
                }

                var usedEnsembl = false
                if (bestMatch == null)
                {
                    if (ensemblGene != null)
                    {
                        sLogger.error("gene: {}*{}, no full match yet has ensembl", geneData.geneName, geneData.allele)

                        if (geneData.region == IgTcrRegion.D_REGION)
                        {
                            // use ensembl for D region
                            geneData.genomicLocationV38 = toGenomicLocation(ensemblGene)
                            usedEnsembl = true
                        }
                    }
                    sLogger.info("gene: {}*{}, no full match", geneData.geneName, geneData.allele)
                }
                else
                {
                    geneData.genomicLocationV38 = matchToQueryGenomicLocation(bestMatch)
                    sLogger.info("gene: {}*{}, match: {}, gene loc: {}", geneData.geneName, geneData.allele, bestMatch, geneData.genomicLocationV38)
                }

                val ciderLocation = if (usedEnsembl) {
                    val chromosomeNum = stripChrPrefix(ensemblGene.Chromosome)
                    "Homo sapiens chromosome ${chromosomeNum}, GRCh38.p13 Primary Assembly:${ensemblGene.GeneStart}-${ensemblGene.GeneEnd}"
                }
                else {
                    if (bestMatch == null) "" else {
                        val start = if (bestMatch.subjectFrame == Strand.FORWARD) bestMatch.subjectAlignStart else bestMatch.subjectAlignEnd
                        val end = if (bestMatch.subjectFrame == Strand.FORWARD) bestMatch.subjectAlignEnd else bestMatch.subjectAlignStart
                        assert(start <= end)
                        "${bestMatch.subjectTitle}:$start-$end"
                    }
                }
                println("GeneBlastnData:${geneData.geneName}\t${geneData.allele}\t${geneData.region}\t${geneData.functionality}\t${geneData.sequenceWithoutGaps}\t$k\t$ciderLocation")
            }
        }

        fun convertGeneTo37(gene: IgTcrGene, genomicLiftOver: LiftOver) : IgTcrGene
        {
            var geneLocationV37: GenomicLocation? = if (gene.geneLocation != null)
                convertGenomicLocationTo37(gene.geneLocation, genomicLiftOver)
            else null

            var anchorLocation: GenomicLocation? = if (gene.anchorLocation != null)
                convertGenomicLocationTo37(gene.anchorLocation, genomicLiftOver)
            else null

            // apply blacklist
            if (gene.geneName in liftOverBlacklist)
            {
                sLogger.info("gene: {} in liftover blacklist, clearing v37 genomic location", gene.geneName)
                geneLocationV37 = null
                anchorLocation = null
            }

            // some anchor lengths are not the same, but we cannot do much about them
            if (anchorLocation != null && anchorLocation.baseLength() != gene.anchorLocation!!.baseLength())
            {
                // decide whether to change anchor start or anchor end
                // the logic is that for V gene we need to preserve the TGT at the end of anchor
                // for J gene we need to preserve the TGG at the start of anchor
                if ((gene.region == IgTcrRegion.V_REGION) == (anchorLocation.strand == Strand.FORWARD))
                {
                    sLogger.warn("gene: {}*{}, strand: {}, different base lengths between v38({}) and v37({}), anchor length: {}, changing anchor start",
                        gene.geneName, gene.allele, anchorLocation.strand, gene.anchorLocation.baseLength(), anchorLocation.baseLength(),
                        gene.anchorSequence!!.length)

                    anchorLocation = anchorLocation.copy(posStart = anchorLocation.posEnd - gene.anchorSequence.length + 1)
                }
                else
                {
                    sLogger.warn("gene: {}*{}, strand: {}, different base lengths between v38({}) and v37({}), anchor length: {}, changing anchor end",
                        gene.geneName, gene.allele, anchorLocation.strand, gene.anchorLocation.baseLength(), anchorLocation.baseLength(),
                        gene.anchorSequence!!.length)

                    anchorLocation = anchorLocation.copy(posEnd = anchorLocation.posStart + gene.anchorSequence.length - 1)
                }
            }

            return gene.copy(geneLocation = geneLocationV37, anchorLocation = anchorLocation)
        }

        fun convertGenomicLocationTo37Hmf(genomicLocV38: GenomicLocation, genomicLiftOver: GenomeLiftoverCache): GenomicLocation?
        {
            if (genomicLocV38.inPrimaryAssembly)
            {
                val v37PosStart: Int = genomicLiftOver.convertPositionTo37(genomicLocV38.chromosome, genomicLocV38.posStart)
                val v37PosEnd: Int = genomicLiftOver.convertPositionTo37(genomicLocV38.chromosome, genomicLocV38.posEnd)

                if (v37PosStart != -1 && v37PosEnd != -1)
                {
                    return if (v37PosEnd < v37PosStart)
                    {
                        // different strand
                        genomicLocV38.copy(chromosome = RefGenomeVersion.V37.versionedChromosome(genomicLocV38.chromosome),
                            posStart = v37PosEnd, posEnd = v37PosStart, strand = genomicLocV38.strand.opposite)
                    }
                    else
                    {
                        genomicLocV38.copy(chromosome = RefGenomeVersion.V37.versionedChromosome(genomicLocV38.chromosome),
                            posStart = v37PosStart, posEnd = v37PosEnd)
                    }
                }
            }
            return null
        }

        // To convert from v38 to v37, we use two methods. First we use the Hartwig GenomicLiftOverCache, and then
        // we use the htsjdk Liftover with downloaded Chain file.
        fun convertGenomicLocationTo37(genomicLocation: GenomicLocation, genomicLiftOver: LiftOver): GenomicLocation?
        {
            if (genomicLocation.inPrimaryAssembly)
            {
                // try the HMF one
                val genomicLocV37 = convertGenomicLocationTo37Hmf(genomicLocation, genomicLiftOverHmf)

                if (genomicLocV37 != null)
                {
                    return genomicLocV37
                }
                else
                {
                    val interval38 = Interval(genomicLocation.chromosome, genomicLocation.posStart, genomicLocation.posEnd, genomicLocation.strand == Strand.REVERSE, "")
                    val interval37 = genomicLiftOver.liftOver(interval38)

                    if (interval37 != null)
                    {
                        sLogger.info("HMF cannot convert but htsjdk liftover can: v38({}) v37({})", genomicLocation, genomicLocV37)
                        return genomicLocation.copy(chromosome = RefGenomeVersion.V37.versionedChromosome(genomicLocation.chromosome),
                            posStart = interval37.start, posEnd = interval37.end, strand = if (interval37.isPositiveStrand) Strand.FORWARD else Strand.REVERSE)
                    }
                }
            }
            return null
        }

        fun toGenomicLocation(ensemblGene: GeneData) : GenomicLocation
        {
            return GenomicLocation(ensemblGene.Chromosome, ensemblGene.GeneStart,
                ensemblGene.GeneEnd, Strand.valueOf(ensemblGene.Strand.toInt()))
        }

        fun matchToQueryGenomicLocation(match: BlastnMatch) : GenomicLocation
        {
            // we need to correct for the ends to make sure things align properly
            val startExtend = match.queryAlignStart - 1
            val endExtend = match.querySeqLen - match.queryAlignEnd
            val matchGenomicLoc = AlignmentUtil.toGenomicLocation(match)!!

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

        fun validateAnchorLocations(igTcrGeneList: List<IgTcrGene>, refGenomeFastaPath: String)
        {
            val refGenomeFile = IndexedFastaSequenceFile(File(refGenomeFastaPath))
            val genomicLocationValidator = GenomicLocationValidator(refGenomeFile)

            for (gene in igTcrGeneList)
            {
                if (gene.anchorLocation != null && gene.anchorLocation.inPrimaryAssembly)
                {
                    if (!genomicLocationValidator.validateAgainstRefGenome(gene.anchorSequence!!, gene.anchorLocation))
                    {
                        sLogger.error("gene: {} anchor location: {} does not match anchor seq: {}",
                            gene.geneAllele, gene.anchorLocation, gene.anchorSequence)
                        // throw RuntimeException()
                    }
                }
            }
        }
    }
}

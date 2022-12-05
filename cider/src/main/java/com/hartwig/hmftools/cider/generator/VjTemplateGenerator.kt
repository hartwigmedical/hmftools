package com.hartwig.hmftools.cider.generator

import com.beust.jcommander.JCommander
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.UnixStyleUsageFormatter
import com.hartwig.hmftools.cider.*
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache
import com.hartwig.hmftools.common.gene.GeneData
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.utils.FileWriterUtils
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import com.hartwig.hmftools.common.utils.config.RefGenomeVersionConverter
import htsjdk.samtools.reference.FastaSequenceFile
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.samtools.util.SequenceUtil.reverseComplement
import org.apache.commons.csv.CSVFormat
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.IOException
import java.util.regex.Pattern
import kotlin.system.exitProcess

const val ANCHOR_DNA_LENGTH: Int = 30
const val KDE_ANCHOR_START_OFFSET: Int = 15
typealias Column = CiderConstants.VjAnchorTemplateTsvColumn

// simple utility to find the anchor sequence for each gene from IMGT fasta data
// example usage:
// java -cp VjTemplateGeneWriter
// -imgt IMGT_genes.fasta
// -ref_genome_version 37
// -ref_genome /data/resources/bucket/reference_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta
// -ensembl_data_dir /data/resources/public/ensembl_data_cache/37
// -kde_region 2:89131735-89132285(-)
// -output igtcr_anchor.37.tsv
class VjTemplateGeneWriter
{
    @Parameter(names = ["-imgt"], required = true, description = "Input IMGT fasta file")
    lateinit var inputImgtFasta: String

    @Parameter(
        names = ["-" + RefGenomeVersion.REF_GENOME_VERSION],
        required = true,
        description = RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC,
        converter = RefGenomeVersionConverter::class)
    lateinit var refGenomeVersion: RefGenomeVersion

    @Parameter(names = ["-ref_genome"], required = true, description = "Reference genome fasta file")
    lateinit var refGenome: String

    @Parameter(names = ["-" + EnsemblDataCache.ENSEMBL_DATA_DIR], required = true, description = EnsemblDataCache.ENSEMBL_DATA_DIR_CFG)
    lateinit var ensemblDataDir: String

    @Parameter(names = ["-kde_region"], required = true, description = "KDE genome region, format is ch:start-end(strand) e.g. 2:89131735-89132285(-)")
    lateinit var kdeRegion: String

    @Parameter(names = ["-output"], required = true, description = "Output TSV file")
    lateinit var outputTsv: String

    fun run(): Int
    {
        val refGenomeFile = IndexedFastaSequenceFile(File(refGenome))

        val imgtFastaFile = FastaSequenceFile(File(inputImgtFasta), false)

        val ensemblDataCache = EnsemblDataCache(ensemblDataDir, refGenomeVersion)
        val ensemblLoadOk = ensemblDataCache.load(true)

        if (!ensemblLoadOk)
        {
            sLogger.error("Ensembl data cache load failed")
            throw RuntimeException("Ensembl data cache load failed")
        }

        val VJAnchorTemplateList: MutableList<VJAnchorTemplate> = ArrayList()

        while (true)
        {
            val sequence = imgtFastaFile.nextSequence()

            if (sequence == null)
                break

            val geneAllele = sequence.name.split('|')[1]
            val toks = geneAllele.split('*')

            if (toks.size < 2)
            {
                sLogger.info("skipping gene: ${geneAllele}")
                continue
            }

            val geneName = toks[0]
            val allele = toks[1]

            val seqStringWithGaps = sequence.baseString!!.uppercase()
            val seqString = seqStringWithGaps.replace(".", "")

            // we only need to worry about V and J genes, so we query the IGHV
            val vjGeneType: VJGeneType

            try
            {
                vjGeneType = if (geneName == "IGKKDE") VJGeneType.IGKKDE else VJGeneType.valueOf(geneName.take(4))
            }
            catch (e: IllegalArgumentException)
            {
                // don't worry about it if it is not a valid type
                continue
            }

            sLogger.info("processing gene: {}", geneAllele)

            // find gene location
            val geneEnsemblData: GeneData? = ensemblDataCache.getGeneDataByName(geneName)
            var geneLoc: GenomeRegionStrand? = null

            if (geneEnsemblData != null)
            {
                geneLoc = GenomeRegionStrand(
                    geneEnsemblData.Chromosome, geneEnsemblData.GeneStart, geneEnsemblData.GeneEnd,
                    if (geneEnsemblData.forwardStrand()) Strand.FORWARD else Strand.REVERSE
                )

                // sLogger.info("gene: ${geneName}, allele: ${allele}, geneLoc: ${geneLoc}")

                // we want to fix up the gene location if there is mismatch
                val refSeq = queryRefSequence(refGenomeFile, geneLoc)
                geneLoc = correctGeneLocation(seqString, refSeq, geneLoc, 20, 30)
            }

            var vjAnchorTemplate: VJAnchorTemplate? = null

            if (vjGeneType.vj == VJ.V)
            {
                // for V gene it should be

                // some sequences are longer
                var anchorIndex: Int = 282
                var anchor = if (seqStringWithGaps.length > anchorIndex) seqStringWithGaps.substring(anchorIndex, Math.min(312, seqStringWithGaps.length)) else ""

                // anchor needs to be at least 27 bases long
                if (anchor.length < 27)
                {
                    // skip this one
                    continue
                }

                // anchor cannot contain .
                if (anchor.contains("."))
                    continue

                // now we find the anchor index again using the one without gap
                anchorIndex = seqString.indexOf(anchor)

                // if the anchor is < 30 bases long, we fill it in with the last C which is TGT
                // it is missing in some
                if (anchor.length < ANCHOR_DNA_LENGTH)
                {
                    anchor = anchor.take(27) + "TGT"
                }

                // for v gene the anchor is near the end, due to introns we much find the sequence from the end backwards
                val anchorOffsetFromEnd = seqString.length - anchorIndex - anchor.length

                var anchorLocation: GenomeRegionStrand? = null

                if (geneLoc != null && anchorIndex != -1)
                {
                    if (geneLoc.strand == Strand.FORWARD)
                    {
                        anchorLocation = GenomeRegionStrand(
                            geneLoc.chromosome,
                            geneLoc.posEnd - anchorOffsetFromEnd - anchor.length + 1,
                            geneLoc.posEnd - anchorOffsetFromEnd, geneLoc.strand)
                    }
                    else
                    {
                        anchorLocation = GenomeRegionStrand(
                            geneLoc.chromosome,
                            geneLoc.posStart + anchorOffsetFromEnd,
                            geneLoc.posStart + anchorOffsetFromEnd + anchor.length - 1, geneLoc.strand)
                    }
                }

                val aaSeq = Codons.aminoAcidFromBases(anchor)

                // v gene
                sLogger.info("IGHV gene: {}, anchor: {}, offset from end: {}, anchor AA: {}", geneName, anchor, anchorOffsetFromEnd, aaSeq)

                vjAnchorTemplate = VJAnchorTemplate(vjGeneType, geneName, allele, geneLoc, seqString, anchor, anchorLocation)
            }

            if (vjGeneType.vj == VJ.J)
            {
                // J gene rules
                // 30 base sequence starting with TGGGG (J-TRP)
                val signatures = if (geneName.startsWith("IGHJ")) listOf("TGGGG") else listOf("TTTG", "TTCG")
                val anchorIndex: Int = seqString.indexOfAny(signatures)
                var anchor = ""

                if (anchorIndex == -1)
                {
                    sLogger.warn("J gene: {}, {} not found, cannot find anchor", geneName, signatures)
                }
                else
                {
                    // some sequences are longer
                    anchor = seqString.substring(anchorIndex).take(ANCHOR_DNA_LENGTH)

                    // anchor needs to be at least 28 bases long
                    //if (anchor.length < 28)

                    // cannot contain .
                    if (anchor.contains("."))
                        continue

                    val aaSeq = Codons.aminoAcidFromBases(anchor)
                    // v gene
                    sLogger.info("IGHJ gene: {}, anchor: {}, anchor AA: {}", geneName, anchor, aaSeq)
                }

                var anchorLocation: GenomeRegionStrand? = null

                if (geneLoc != null && anchorIndex != -1)
                {
                    if (geneLoc.strand == Strand.FORWARD)
                    {
                        anchorLocation = GenomeRegionStrand(
                            geneLoc.chromosome,
                            geneLoc.posStart + anchorIndex,
                            geneLoc.posStart + anchorIndex + anchor.length - 1, geneLoc.strand)
                    }
                    else
                    {
                        anchorLocation = GenomeRegionStrand(
                            geneLoc.chromosome,
                            geneLoc.posEnd - anchorIndex - anchor.length + 1,
                            geneLoc.posEnd - anchorIndex, geneLoc.strand)
                    }
                }

                vjAnchorTemplate = VJAnchorTemplate(vjGeneType, geneName, allele, geneLoc, seqString, anchor, anchorLocation)
            }

            if (vjAnchorTemplate != null)
            {
                if (vjAnchorTemplate.anchorLocation != null)
                {
                    // validate the anchor sequence
                    validateAgainstRefGenome(vjAnchorTemplate.anchorSequence, vjAnchorTemplate.anchorLocation!!, refGenomeFile)
                }
                VJAnchorTemplateList.add(vjAnchorTemplate)
            }
        }

        // also add the KDE
        VJAnchorTemplateList.add(createKdeRegionTemplate(kdeRegion, refGenomeFile))

        writeOutput(outputTsv, VJAnchorTemplateList)

        return 0
    }

    companion object
    {
        private val sLogger = LogManager.getLogger(VjTemplateGeneWriter::class.java)

        @JvmStatic
        fun main(args: Array<String>)
        {
            // here we have some voodoo to work out if we are being used in pipeline mode or the standalone mode
            val app = VjTemplateGeneWriter()
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

        @Throws(IOException::class)
        fun writeOutput(filename: String, VJAnchorTemplateList: List<VJAnchorTemplate>)
        {
            val csvFormat = CSVFormat.Builder.create()
                .setDelimiter('\t').setRecordSeparator('\n')
                .setHeader(Column::class.java)
                .build()

            csvFormat.print(FileWriterUtils.createBufferedWriter(filename)).use { csvPrinter ->
                for (gene: VJAnchorTemplate in VJAnchorTemplateList)
                {
                    for (col in Column.values())
                    {
                        when (col)
                        {
                            Column.gene -> csvPrinter.print(gene.geneName)
                            Column.allele -> csvPrinter.print(gene.allele)
                            Column.chr -> csvPrinter.print(if (gene.geneLocation != null) gene.geneLocation.chromosome else "")
                            Column.posStart -> csvPrinter.print(if (gene.geneLocation != null) gene.geneLocation.posStart else -1)
                            Column.posEnd -> csvPrinter.print(if (gene.geneLocation != null) gene.geneLocation.posEnd else -1)
                            Column.strand -> csvPrinter.print(if (gene.geneLocation != null) gene.geneLocation.strand.asChar() else "")
                            Column.anchorStart -> csvPrinter.print(if (gene.anchorLocation != null) gene.anchorLocation.posStart else -1)
                            Column.anchorEnd -> csvPrinter.print(if (gene.anchorLocation != null) gene.anchorLocation.posEnd else -1)
                            Column.anchorSequence -> csvPrinter.print(gene.anchorSequence)
                            Column.anchorAA -> csvPrinter.print(Codons.aminoAcidFromBases(gene.anchorSequence))
                            Column.sequence -> csvPrinter.print(gene.sequence)
                        }
                    }
                    csvPrinter.println()
                }
            }
        }

        fun correctGeneLocation(seq: String, refSeq: String,
                                refGenomeRegionStrand: GenomeRegionStrand, minBasesToCompare: Int, maxBasesToCompare: Int)
            : GenomeRegionStrand
        {
            // the sequence from the two files are not exact matches, so we need to correct for it
            // note that we only try to find the first section of the sequence in ref seq
            // this is to ensure it works even if there are some bases inserted / deleted somewhere
            val (startShift, numStartMismatch) = calcPositionMismatch(seq, refSeq, minBasesToCompare, maxBasesToCompare)
            var (endShift, numEndMismatch) = calcPositionMismatch(seq.reversed(), refSeq.reversed(), minBasesToCompare, maxBasesToCompare)
            endShift = -endShift

            if (numStartMismatch >= 10 || numEndMismatch >= 10)
            {
                // the sequences between ref and the location we got from IMGT is not a full match. This would be ok as long as
                // the anchor location is a good match
                sLogger.warn("mismatch count({}) non matching sequence: seq({}) refSeq({})", numStartMismatch, seq, refSeq)
            }

            if (startShift == 0 && seq.length == refSeq.length)
                return refGenomeRegionStrand

            val correctedGeneLoc = if (refGenomeRegionStrand.strand == Strand.FORWARD)
                GenomeRegionStrand(refGenomeRegionStrand.chromosome, refGenomeRegionStrand.posStart + startShift,
                    refGenomeRegionStrand.posEnd + endShift, refGenomeRegionStrand.strand)
            else
                GenomeRegionStrand(refGenomeRegionStrand.chromosome, refGenomeRegionStrand.posStart - endShift,
                    refGenomeRegionStrand.posEnd - startShift, refGenomeRegionStrand.strand)

            sLogger.info("seq({}) refSeq({}) refLocation({}) shift({}, {}) corrected({})",
                seq, refSeq, refGenomeRegionStrand, startShift, endShift, correctedGeneLoc)

            return correctedGeneLoc
        }

        fun calcPositionMismatch(seq: String, refSeq: String, minBasesToCompare: Int, maxBasesToCompare: Int): Pair<Int, Int>
        {
            val minBasesToCompare = minOf(seq.length, refSeq.length, minBasesToCompare)

            var bestShift = 0
            var lowestNumMismatch = Int.MAX_VALUE

            val minShift = -seq.length + minBasesToCompare
            val maxShift = refSeq.length - minBasesToCompare

            // we shift along to find the best match
            for (i in minShift .. maxShift)
            {
                var j = Math.max(0, i)
                var numMismatch = 0
                var numCompare = 0
                while ((j - i) < seq.length && j < refSeq.length && numCompare < maxBasesToCompare)
                {
                    if (seq[j - i] != refSeq[j])
                    {
                        ++numMismatch
                    }
                    ++numCompare
                    ++j
                }

                if (numCompare >= minBasesToCompare && numMismatch < lowestNumMismatch)
                {
                    lowestNumMismatch = numMismatch
                    bestShift = i
                }
            }

            return Pair(bestShift, lowestNumMismatch)
        }

        // validate sequence against the ref genome file to make sure we got it right
        fun validateAgainstRefGenome(seq: String, genomeRegionStrand: GenomeRegionStrand, refGenome: IndexedFastaSequenceFile) : Boolean
        {
            val refGenomeSeq = queryRefSequence(refGenome, genomeRegionStrand)

            if (refGenomeSeq.length != seq.length)
            {
                sLogger.warn("validation failed: seq({}) and ref genome seq({} of {}) length mismatch", seq, refGenomeSeq, genomeRegionStrand)
                return false
            }

            // do not allow more than 2 bases difference
            var numDiff: Int = 0
            for (i in seq.indices)
            {
                if (seq[i] != refGenomeSeq[i])
                    ++numDiff
            }

            if (numDiff > 6)
            {
                sLogger.error("validation failed: seq({}) and ref genome seq({} of {}) sequence mismatch({}) > 6", seq, refGenomeSeq, genomeRegionStrand, numDiff)
                return false
            }
            return true
        }

        fun queryRefSequence(refGenome: IndexedFastaSequenceFile, genomeRegionStrand: GenomeRegionStrand): String
        {
            var chromosome = genomeRegionStrand.chromosome

            if (!refGenome.index.hasIndexEntry(chromosome))
            {
                // maybe need to try removing chr
                chromosome = chromosome.replace("chr", "")
            }

            val forwardSeq = refGenome.getSubsequenceAt(chromosome,
                genomeRegionStrand.start().toLong(), genomeRegionStrand.end().toLong()).baseString
            if (genomeRegionStrand.strand == Strand.FORWARD)
                return forwardSeq
            else
                return reverseComplement(forwardSeq)
        }

        fun createKdeRegionTemplate(kdeRegion: String, refGenome: IndexedFastaSequenceFile) : VJAnchorTemplate
        {
            // format of kde region string is 2:89131735-89132285(-)
            val genomeRegionPattern = Pattern.compile("""(.+):(\d+)-(\d+)\((.)\)""")
            val matcher = genomeRegionPattern.matcher(kdeRegion)
            if (!matcher.matches())
            {
                sLogger.error("input KDE region string: {} parse failed", kdeRegion)
                throw RuntimeException("input KDE region string: ${kdeRegion} parse failed")
            }

            val chr = matcher.group(1)
            val start = matcher.group(2).toInt()
            val end = matcher.group(3).toInt()
            val strand = matcher.group(4)[0]

            val genomeRegionStrand = GenomeRegionStrand(chr, start, end, Strand.valueOf(strand))
            val kdeSeq = queryRefSequence(refGenome, genomeRegionStrand)

            // get some sort of "anchor", in reality does not really matter
            val anchorGenomeRegionStrand =
            // now the anchor genome region, depends on the strand
            if (genomeRegionStrand.strand == Strand.FORWARD)
            {
                val posStart = genomeRegionStrand.posStart + KDE_ANCHOR_START_OFFSET
                genomeRegionStrand.copy(
                    posStart = posStart,
                    posEnd = posStart + ANCHOR_DNA_LENGTH - 1
                )
            }
            else
            {
                val posEnd = genomeRegionStrand.posEnd - KDE_ANCHOR_START_OFFSET
                genomeRegionStrand.copy(
                    posStart = posEnd - ANCHOR_DNA_LENGTH + 1,
                    posEnd = posEnd
                )
            }

            val anchorSeq = queryRefSequence(refGenome, anchorGenomeRegionStrand)

            // now we can create the template
            return VJAnchorTemplate(VJGeneType.IGKKDE, "IGKKDE", "01",
                genomeRegionStrand, kdeSeq, anchorSeq, anchorGenomeRegionStrand)
        }
    }
}

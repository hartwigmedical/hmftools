package com.hartwig.hmftools.cider

import com.beust.jcommander.JCommander
import com.beust.jcommander.Parameter
import com.beust.jcommander.ParameterException
import com.beust.jcommander.UnixStyleUsageFormatter
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.genome.region.Strand
import com.hartwig.hmftools.common.sequence.SequenceAligner
import com.hartwig.hmftools.common.utils.FileWriterUtils
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator
import htsjdk.samtools.reference.FastaSequenceFile
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.samtools.util.SequenceUtil.reverseComplement
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager
import java.io.*
import kotlin.system.exitProcess

// simple utility to find the anchor sequence for each gene
// example usage:
// -imgt human_IMGT+C.fa
// -gene_loc hg19_bcrtcr.fa
// -ref_genome /data/resources/bucket/reference_genome/37/Homo_sapiens.GRCh37.GATK.illumina.fasta
// -output jv_anchor.tsv
class VjTemplateGeneWriter
{
    @Parameter(names = ["-imgt"], required = true, description = "Input IMGT fasta file")
    var inputImgtFasta: String = ""

    @Parameter(names = ["-gene_loc"], required = true, description = "Input gene location fasta file")
    var inputGeneLocFasta: String = ""

    @Parameter(names = ["-ref_genome"], required = true, description = "Reference genome fasta file")
    var refGenome: String = ""

    @Parameter(names = ["-output"], required = true, description = "Output TSV file")
    var outputTsv: String = ""

    fun run(): Int
    {
        val refGenomeFile = IndexedFastaSequenceFile(File(refGenome))

        val imgtFaFile = FastaSequenceFile(File(inputImgtFasta), false)

        val geneLocMap: Map<String, Pair<GeneLocation, String>> = geneLocationMap()

        val VJGeneList: MutableList<VJGene> = ArrayList()

        while (true)
        {
            val sequence = imgtFaFile.nextSequence()

            if (sequence == null)
                break

            //sLogger.info(sequence.name)

            val toks = sequence.name.split('*')

            if (toks.size < 2)
            {
                sLogger.info("skipping gene: ${sequence.name}")
                continue
            }

            val geneName = toks[0]
            val allele = toks[1]

            // find gene location
            var (geneLoc: GeneLocation?, geneLocSeq: String?) = geneLocMap.getOrDefault(geneName, Pair(null, null))

            // sLogger.info("gene: ${geneName}, allele: ${allele}, geneLoc: ${geneLoc}")

            val seqStringWithGaps = sequence.baseString!!
            val seqString = seqStringWithGaps.replace(".", "")

            // we want to fix up the gene location if there is mismatch
            if (geneLoc != null && geneLocSeq != null)
            {
                geneLoc = correctGeneLocation(seqString, geneLocSeq, geneLoc, 20, 30)
            }

            var VJGene: VJGene? = null

            // if (geneLoc != null)

            if (geneName.startsWith("IGHV") ||
                geneName.startsWith("IGKV") ||
                geneName.startsWith("IGLV") ||
                geneName.startsWith("TRAV") ||
                geneName.startsWith("TRBV") ||
                geneName.startsWith("TRDV") ||
                geneName.startsWith("TRGV"))
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
                if (anchor.length < 30)
                {
                    anchor = anchor.take(27) + "TGT"
                }

                // for v gene the anchor is near the end, due to introns we much find the sequence from the end backwards
                val anchorOffsetFromEnd = seqString.length - anchorIndex - anchor.length

                var anchorLocation: GeneLocation? = null

                if (geneLoc != null && anchorIndex != -1)
                {
                    if (geneLoc.strand == Strand.FORWARD)
                    {
                        anchorLocation = GeneLocation(
                            geneLoc.chromosome,
                            geneLoc.posEnd - anchorOffsetFromEnd - anchor.length + 1,
                            geneLoc.posEnd - anchorOffsetFromEnd, geneLoc.strand)
                    }
                    else
                    {
                        anchorLocation = GeneLocation(
                            geneLoc.chromosome,
                            geneLoc.posStart + anchorOffsetFromEnd,
                            geneLoc.posStart + anchorOffsetFromEnd + anchor.length - 1, geneLoc.strand)
                    }
                }

                val aaSeq = Codons.aminoAcidFromBases(anchor)

                // v gene
                sLogger.info("IGHV gene: {}, anchor: {}, offset from end: {}, anchor AA: {}", geneName, anchor, anchorOffsetFromEnd, aaSeq)

                VJGene = VJGene(sequence.name, geneName, allele, geneLoc, seqString, anchor, anchorLocation)
            }

            if (geneName.startsWith("IGHJ") ||
                geneName.startsWith("IGKJ") ||
                geneName.startsWith("IGLJ") ||
                geneName.startsWith("TRAJ") ||
                geneName.startsWith("TRBJ") ||
                geneName.startsWith("TRDJ") ||
                geneName.startsWith("TRGJ"))
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
                    anchor = seqString.substring(anchorIndex).take(30)

                    // anchor needs to be at least 28 bases long
                    //if (anchor.length < 28)

                    // cannot contain .
                    if (anchor.contains("."))
                        continue

                    val aaSeq = Codons.aminoAcidFromBases(anchor)
                    // v gene
                    sLogger.info("IGHJ gene: {}, anchor: {}, anchor AA: {}", geneName, anchor, aaSeq)
                }

                var anchorLocation: GeneLocation? = null

                if (geneLoc != null && anchorIndex != -1)
                {
                    if (geneLoc.strand == Strand.FORWARD)
                    {
                        anchorLocation = GeneLocation(
                            geneLoc.chromosome,
                            geneLoc.posStart + anchorIndex,
                            geneLoc.posStart + anchorIndex + anchor.length - 1, geneLoc.strand)
                    }
                    else
                    {
                        anchorLocation = GeneLocation(
                            geneLoc.chromosome,
                            geneLoc.posEnd - anchorIndex - anchor.length + 1,
                            geneLoc.posEnd - anchorIndex, geneLoc.strand)
                    }
                }

                VJGene = VJGene(sequence.name, geneName, allele, geneLoc, seqString, anchor, anchorLocation)
            }

            if (VJGene != null)
            {
                if (VJGene.anchorLocation != null)
                {
                    // validate the anchor sequence
                    validateAgainstRefGenome(VJGene.anchorSequence, VJGene.anchorLocation!!, refGenomeFile)
                }
                VJGeneList.add(VJGene)
            }
        }

        writeOutput(outputTsv, VJGeneList)

        return 0
    }

    // return a map of gene name to their location and the sequence as indicated in the location file
    fun geneLocationMap() : Map<String, Pair<GeneLocation, String>>
    {
        FastaSequenceFile(File(inputGeneLocFasta), false).use { geneLocFaFile ->
            val geneLocMap: MutableMap<String, Pair<GeneLocation, String>> = HashMap()

            while (true)
            {
                val sequence = geneLocFaFile.nextSequence()

                if (sequence == null)
                    break

                // sLogger.info(sequence.name)

                val toks = sequence.name.split(' ')
                val geneName = toks[0]
                val chromosome = toks[1]
                val startPos = toks[2].toInt()
                val endPos = toks[3].toInt()
                val strand = if (toks[4] == "+") Strand.FORWARD else Strand.REVERSE

                geneLocMap[geneName] = Pair(GeneLocation(chromosome, startPos, endPos, strand), sequence.baseString)
            }

            return geneLocMap
        }
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

        enum class Column {
            id,
            gene,
            allele,
            chr,
            posStart,
            posEnd,
            strand,
            anchorStart,
            anchorEnd,
            anchorSequence,
            anchorAA,
            sequence
        }

        @Throws(IOException::class)
        fun writeOutput(filename: String, VJGeneList: List<VJGene>)
        {
            val csvFormat = CSVFormat.Builder.create()
                .setDelimiter('\t').setRecordSeparator('\n')
                .setHeader(Column::class.java)
                .build()

            csvFormat.print(FileWriterUtils.createBufferedWriter(filename)).use { csvPrinter ->
                for (gene: VJGene in VJGeneList)
                {
                    for (col in Column.values())
                    {
                        when (col)
                        {
                            Column.id -> csvPrinter.print(gene.id)
                            Column.gene -> csvPrinter.print(gene.name)
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

        fun correctGeneLocation(seq: String, refSeq: String, refGeneLocation: GeneLocation, minBasesToCompare: Int, maxBasesToCompare: Int) : GeneLocation
        {
            // the sequence from the two files are not exact matches, so we need to correct for it
            // note that we only try to find the first section of the sequence in ref seq
            // this is to ensure it works even if there are some bases inserted / deleted somewhere
            val (startShift, numStartMismatch) = calcPositionMismatch(seq, refSeq, minBasesToCompare, maxBasesToCompare)
            var (endShift, numEndMismatch) = calcPositionMismatch(seq.reversed(), refSeq.reversed(), minBasesToCompare, maxBasesToCompare)
            endShift = -endShift

            if (numStartMismatch >= 10 || numEndMismatch >= 10)
            {
                sLogger.warn("non matching sequence: seq({}) refSeq({}) mismatch count({})", seq, refSeq, numStartMismatch)
            }

            if (startShift == 0 && seq.length == refSeq.length)
                return refGeneLocation

            val correctedGeneLoc = if (refGeneLocation.strand == Strand.FORWARD)
                GeneLocation(refGeneLocation.chromosome, refGeneLocation.posStart + startShift,
                    refGeneLocation.posEnd + endShift, refGeneLocation.strand)
            else
                GeneLocation(refGeneLocation.chromosome, refGeneLocation.posStart - endShift,
                    refGeneLocation.posEnd - startShift, refGeneLocation.strand)

            sLogger.info("seq({}) refSeq({}) refLocation({}) shift({}, {}) corrected({})",
                seq, refSeq, refGeneLocation, startShift, endShift, correctedGeneLoc)

            return correctedGeneLoc
        }

        fun calcPositionMismatch(seq: String, refSeq: String, minBasesToCompare: Int, maxBasesToCompare: Int): Pair<Int, Int>
        {
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
        fun validateAgainstRefGenome(seq: String, geneLocation: GeneLocation, refGenome: IndexedFastaSequenceFile) : Boolean
        {
            val refGenomeSeq = refSequence(refGenome, geneLocation)

            if (refGenomeSeq.length != seq.length)
            {
                sLogger.warn("validation failed: seq({}) and ref genome seq({} of {}) length mismatch", seq, refGenomeSeq, geneLocation)
                return false
            }

            // do not allow more than 2 bases difference
            var numDiff: Int = 0
            for (i in seq.indices)
            {
                if (seq[i] != refGenomeSeq[i])
                    ++numDiff
            }

            if (numDiff >= 5)
            {
                sLogger.error("validation failed: seq({}) and ref genome seq({} of {}) sequence mismatch({}) >= 5", seq, refGenomeSeq, geneLocation, numDiff)
                return false
            }
            return true
        }

        fun refSequence(refGenome: IndexedFastaSequenceFile, geneLocation: GeneLocation): String
        {
            var chromosome = geneLocation.chromosome

            if (!refGenome.index.hasIndexEntry(chromosome))
            {
                // maybe need to try removing chr
                chromosome = chromosome.replace("chr", "");
            }

            val forwardSeq = refGenome.getSubsequenceAt(chromosome,
                geneLocation.start().toLong(), geneLocation.end().toLong()).baseString
            if (geneLocation.strand == Strand.FORWARD)
                return forwardSeq
            else
                return reverseComplement(forwardSeq)
        }

        fun findAnchorPos(row: CSVRecord, refGenome: IndexedFastaSequenceFile): Pair<Long, Long>
        {
            if (row["Chr"].isEmpty())
                return Pair(-1, -1)

            // add error correction term to pos start / pos end just in case
            val ERROR_CORRECT: Long = 200
            val chr: String = row["Chr"].substring(3)
            val posStart: Long = row["Pos Start"].toLong() - ERROR_CORRECT
            val posEnd: Long = row["Pos End"].toLong() + ERROR_CORRECT
            val strand: String = row["Strand"]

            val refSeq = refSequence(refGenome, chr, posStart, posEnd, strand)
            val anchorSeq = row["AnchorSequence (30 bases ending with conserved C)"]

            if (anchorSeq.isEmpty())
                return Pair(-1, -1)

            // now try to use aligner to find it in side
            val alignment: SequenceAligner.Alignment = SequenceAligner.alignSubsequence(anchorSeq, refSeq)

            val alignStart = alignment.alignOps.indexOfFirst({ op -> op == SequenceAligner.AlignOp.MATCH || op == SequenceAligner.AlignOp.SUBSTITUTION})
            val alignEnd = alignment.alignOps.indexOfFirst({ op -> op == SequenceAligner.AlignOp.MATCH || op == SequenceAligner.AlignOp.SUBSTITUTION})

            if (alignStart == -1 || alignEnd == -1)
            {
                sLogger.warn("cannot align: ${row}")
                return Pair(-1, -1)
            }

            if (alignEnd - alignStart + 1 > anchorSeq.length + 4)
            {
                sLogger.warn("cannot align: ${row}, too many insertions")
                return Pair(-1, -1)
            }

            val anchorStart =
            if (row["Strand"] == "-")
                posEnd - alignStart - anchorSeq.length + 1
            else
                posStart + alignStart

            // since we are one based
            val anchorEnd = anchorStart + anchorSeq.length - 1

            sLogger.info("pos = ${alignStart}, anchorStart = ${anchorStart}, anchorEnd = ${anchorEnd}")

            if (row["Gene"] == "IGHV3-32")
            {
                val checkSeq = refGenome.getSubsequenceAt(row["Chr"].substring(3), anchorStart, anchorEnd).baseString
                sLogger.info("here: ${anchorSeq} vs ${checkSeq}, ref seq: ${refSeq}")
            }

            return Pair(anchorStart, anchorEnd)
        }

        fun refSequence(refGenome: IndexedFastaSequenceFile, chr: String, start: Long, end: Long, strand: String): String
        {
            val forwardSeq = refGenome.getSubsequenceAt(chr,
                start, end).baseString
            if (strand == "-")
                return reverseComplement(forwardSeq)
            else
                return forwardSeq
        }
    }
}
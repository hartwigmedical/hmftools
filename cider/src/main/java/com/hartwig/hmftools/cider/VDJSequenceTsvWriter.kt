package com.hartwig.hmftools.cider

import com.hartwig.hmftools.cider.CiderConstants.MIN_NON_SPLIT_READ_STRADDLE_LENGTH
import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.FileWriterUtils
import com.hartwig.hmftools.common.utils.IntPair
import htsjdk.samtools.SAMRecord
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashMap

object VDJSequenceTsvWriter
{
    const val CDR3_FILTER_AA_MIN_LENGTH = 4
    const val CDR3_FILTER_AA_MAX_LENGTH = 40

    private enum class Column
    {
        cdr3Seq,
        cdr3AA,
        filter,
        minHighQualBaseReads,
        assignedReads,
        vAlignedReads,
        jAlignedReads,
        inFrame,
        containsStop,
        vType,
        vAnchorEnd,
        vAnchorSeq,
        vAnchorTemplateSeq,
        vAnchorAA,
        vAnchorTemplateAA,
        vMatchMethod,
        vSimilarityScore,
        vNonSplitReads,
        jType,
        jAnchorStart,
        jAnchorSeq,
        jAnchorTemplateSeq,
        jAnchorAA,
        jAnchorTemplateAA,
        jMatchMethod,
        jSimilarityScore,
        jNonSplitReads,
        layoutId,
        vdjSeq,
        support
    }

    private const val FILE_EXTENSION = ".cider.vdj_seq.tsv"
    private val sLogger = LogManager.getLogger(VDJSequenceTsvWriter::class.java)

    @JvmStatic
    fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writeVDJSequences(
        basePath: String, sample: String, vdjSequences: List<VDJSequence>,
        adaptor: VJReadLayoutAdaptor, reportPartialSeq: Boolean)
    {
        val filePath = generateFilename(basePath, sample)

        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .build()

        var sortedVdj = vdjSequences.sortedWith(
            Collections.reverseOrder(
                Comparator.comparingInt({ vdj: VDJSequence -> vdj.cdr3SupportMin })
                    .thenComparingInt({ vdj: VDJSequence -> vdj.numReads }) // handle the highest quality ones first
                    .thenComparingInt({ vdj: VDJSequence -> vdj.length })
                    .thenComparing({ vdj: VDJSequence -> vdj.cdr3Sequence })
            ))

        if (!reportPartialSeq)
        {
            // remove all the partially rearranged sequences
            sortedVdj = sortedVdj.filter({ seq: VDJSequence -> seq.isFullyRearranged })
        }

        // we create a set of all VDJs that are NOT duplicates
        val notDuplicateVdjs = HashMap<String, VDJSequence>()

        for (vdjSeq in sortedVdj)
        {
            val cdr3 = vdjSeq.cdr3Sequence

            val vdjWithMoreSupport : VDJSequence? = notDuplicateVdjs.get(cdr3)

            if (vdjWithMoreSupport == null)
            {
                notDuplicateVdjs[cdr3] = vdjSeq
            }
            else
            {
                // just make sure it has more support
                require(vdjWithMoreSupport.cdr3SupportMin >= vdjSeq.cdr3SupportMin)
            }
        }

        csvFormat.print(FileWriterUtils.createBufferedWriter(filePath)).use { printer: CSVPrinter ->
            for (vdj in sortedVdj)
            {
                val vdjWithMoreSupport : VDJSequence? = notDuplicateVdjs.get(vdj.cdr3Sequence)
                assert(vdjWithMoreSupport != null)
                val isDuplicate: Boolean = vdjWithMoreSupport !== vdj
                writeVDJSequence(printer, vdj, isDuplicate, adaptor)
            }
        }
    }

    private fun writeVDJSequence(csvPrinter: CSVPrinter, vdj: VDJSequence, isDuplicate: Boolean, adaptor: VJReadLayoutAdaptor)
    {
        val vAnchorByReadMatch: VJAnchorByReadMatch? = vdj.vAnchor as? VJAnchorByReadMatch
        val jAnchorByReadMatch: VJAnchorByReadMatch? = vdj.jAnchor as? VJAnchorByReadMatch

        for (c in Column.values())
        {
            when (c)
            {
                Column.cdr3Seq -> csvPrinter.print(vdj.cdr3Sequence)
                Column.cdr3AA -> csvPrinter.print(cdr3AminoAcidFs(vdj))
                Column.filter -> csvPrinter.print(filterString(vdj, isDuplicate))
                Column.minHighQualBaseReads -> csvPrinter.print(vdj.cdr3SupportMin)
                Column.assignedReads -> csvPrinter.print(vdj.numReads)
                Column.vAlignedReads -> csvPrinter.print(vAnchorByReadMatch?.numReads ?: 0)
                Column.jAlignedReads -> csvPrinter.print(jAnchorByReadMatch?.numReads ?: 0)
                Column.inFrame -> csvPrinter.print(vdj.isInFrame)
                Column.containsStop -> csvPrinter.print(vdj.aminoAcidSequence.contains(Codons.STOP_AMINO_ACID))
                Column.vType -> csvPrinter.print(vdj.vAnchor?.geneType)
                Column.vAnchorEnd -> csvPrinter.print(vdj.vAnchor?.anchorBoundary)
                Column.vAnchorSeq -> csvPrinter.print(vdj.vAnchorSequence)
                Column.vAnchorTemplateSeq -> csvPrinter.print(vdj.vAnchor?.templateAnchorSeq)
                Column.vAnchorAA -> csvPrinter.print(aminoAcidFromBases(vdj.vAnchorSequence))
                Column.vAnchorTemplateAA -> csvPrinter.print(if (vdj.vAnchor != null) aminoAcidFromBases(vdj.vAnchor.templateAnchorSeq) else null)
                Column.vMatchMethod -> csvPrinter.print(vdj.vAnchor?.matchMethod)
                Column.vSimilarityScore -> csvPrinter.print(if (vdj.vAnchor != null) calcAnchorSimilarity(vdj, vdj.vAnchor) else null)
                Column.vNonSplitReads -> csvPrinter.print(countNonSplitReads(vdj, VJ.V, adaptor))
                Column.jType -> csvPrinter.print(vdj.jAnchor?.geneType)
                Column.jAnchorStart -> csvPrinter.print(vdj.jAnchor?.anchorBoundary)
                Column.jAnchorSeq -> csvPrinter.print(vdj.jAnchorSequence)
                Column.jAnchorTemplateSeq -> csvPrinter.print(vdj.jAnchor?.templateAnchorSeq)
                Column.jAnchorAA -> csvPrinter.print(aminoAcidFromBases(vdj.jAnchorSequence))
                Column.jAnchorTemplateAA -> csvPrinter.print(if (vdj.jAnchor != null) aminoAcidFromBases(vdj.jAnchor.templateAnchorSeq) else null)
                Column.jMatchMethod -> csvPrinter.print(vdj.jAnchor?.matchMethod)
                Column.jSimilarityScore -> csvPrinter.print(if (vdj.jAnchor != null) calcAnchorSimilarity(vdj, vdj.jAnchor) else null)
                Column.jNonSplitReads -> csvPrinter.print(countNonSplitReads(vdj, VJ.J, adaptor))
                Column.layoutId -> csvPrinter.print(vdj.layout.id)
                Column.vdjSeq -> csvPrinter.print(vdj.sequence)
                Column.support -> csvPrinter.print(CiderUtils.countsToString(vdj.supportCounts))
            }
        }
        csvPrinter.println()
    }

    // we want to replace X with _, easier to see in file
    fun aminoAcidFromBases(dna: String): String
    {
        return Codons.aminoAcidFromBases(dna).replace(Codons.STOP_AMINO_ACID, '_')
    }

    // add suffix for out of frame
    fun cdr3AminoAcidFs(vdj: VDJSequence): String
    {
        val suffix = if (vdj.isInFrame) "" else "fs"
        return aminoAcidFromBases(vdj.cdr3Sequence) + suffix
    }

    fun calcAnchorSimilarity(vdj: VDJSequence, anchor: VJAnchor) : Int
    {
        val seq: String
        val templateAnchorSeq: String = anchor.templateAnchorSeq

        if (anchor.vj == VJ.V)
        {
            seq = vdj.vAnchorSequence
        }
        else
        {
            seq = vdj.jAnchorSequence
        }

        try
        {
            return BlosumSimilarityCalc.calcSimilarityScore(anchor.vj, templateAnchorSeq, seq)
        }
        catch (e: java.lang.IllegalArgumentException)
        {
            throw IllegalArgumentException("cannot calc similarity score: ${templateAnchorSeq} and ${seq}")
        }
    }

    fun filterString(vdj: VDJSequence, isDuplicate: Boolean): String
    {
        val filters = ArrayList<String>()

        if (vdj.vAnchor == null)
        {
            filters.add("NO_V_ANCHOR")
        }
        if (vdj.jAnchor == null)
        {
            filters.add("NO_J_ANCHOR")
        }
        if ((vdj.vAnchor is VJAnchorByBlosum) && (vdj.vAnchor.similarityScore < 0))
        {
            filters.add("POOR_V_ANCHOR")
        }
        if ((vdj.jAnchor is VJAnchorByBlosum) && (vdj.jAnchor.similarityScore < 0))
        {
            filters.add("POOR_J_ANCHOR")
        }
        if (isDuplicate)
        {
            filters.add("DUPLICATE")
        }
        if (vdj.cdr3Sequence.length < CDR3_FILTER_AA_MIN_LENGTH * 3)
        {
            filters.add("MIN_LENGTH")
        }
        if (vdj.cdr3Sequence.length > CDR3_FILTER_AA_MAX_LENGTH * 3)
        {
            filters.add("MAX_LENGTH")
        }
        if (filters.isEmpty())
            filters.add("PASS")

        return filters.joinToString(separator = ";")
    }

    // for now we want to just log where the mappings are
    fun countNonSplitReads(vdj: VDJSequence, vj: VJ, adaptor: VJReadLayoutAdaptor) : Int
    {
        val alignedPos = vdj.layout.alignedPosition - vdj.layoutSliceStart

        if (vj == VJ.V && vdj.vAnchor == null)
            return 0

        if (vj == VJ.J && vdj.jAnchor == null)
            return 0

        val boundaryPos: Int = if (vj == VJ.V)
        {
            vdj.vAnchorBoundary!!
        }
        else
        {
            vdj.jAnchorBoundary!!
        }
        var nonSplitReadCount = 0

        for (read in vdj.layout.reads)
        {
            val samRecord: SAMRecord = adaptor.toReadCandidate(read).read
            val layoutReadSlice: ReadSlice = adaptor.toLayoutReadSlice(read)

            // AGATCTGAG-GACACGGCCGTGTATTACTGT-GCGAGAGACACAGTGTGAAAACCCACATCCTGAGAGTGTCAGAAACCCTGAGGGA
            //           |___________________|
            //                 V anchor
            //                  =========================> read slice
            //                                        |  <-- aligned position
            //                  |------------|           <-- this is the value we want
            val readPosWithinVdj = alignedPos - read.alignedPosition
            val readSliceAnchorBoundary = boundaryPos - readPosWithinVdj

            // work out where in the VDJ sequence is this read mapped
            if (!samRecord.readUnmappedFlag)
            {
                for (alignBlock in CiderUtils.getAdjustedAlignmentBlocks(samRecord.cigar))
                {
                    // now get those positions in terms of read slice
                    val alignRangeInReadSlice: IntPair = layoutReadSlice.readRangeToSliceRange(
                        alignBlock.readStart - 1,
                        alignBlock.readStart - 1 + alignBlock.length)

                    require(alignRangeInReadSlice.left < alignRangeInReadSlice.right)

                    // to count as a non split read, it needs to be away from the boundary
                    if (readSliceAnchorBoundary >= alignRangeInReadSlice.left + MIN_NON_SPLIT_READ_STRADDLE_LENGTH &&
                        readSliceAnchorBoundary <= alignRangeInReadSlice.right - MIN_NON_SPLIT_READ_STRADDLE_LENGTH)
                    {
                        ++nonSplitReadCount

                        // this read straddles a v anchor boundary
                        sLogger.debug("read({}) cigar({}) revcomp({}), straddles {} boundary, align offset({}:{}), boundary offset({})",
                            samRecord, samRecord.cigarString, layoutReadSlice.reverseComplement, vj,
                            alignRangeInReadSlice.left, alignRangeInReadSlice.right, readSliceAnchorBoundary)
                    }
                }
            }
        }
        return nonSplitReadCount
    }
}
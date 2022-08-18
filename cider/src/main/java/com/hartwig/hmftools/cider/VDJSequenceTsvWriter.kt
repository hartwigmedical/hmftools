package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
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
        vSimilarityScore,
        jType,
        jAnchorStart,
        jAnchorSeq,
        jAnchorTemplateSeq,
        jAnchorAA,
        jAnchorTemplateAA,
        jSimilarityScore,
        layoutId,
        vdjSeq,
        support
    }

    private const val FILE_EXTENSION = ".cider.vdj_seq.tsv"

    @JvmStatic
    fun generateFilename(basePath: String, sample: String): String
    {
        return basePath + File.separator + sample + FILE_EXTENSION
    }

    @JvmStatic
    fun writeVDJSequences(basePath: String, sample: String, vdjSequences: List<VDJSequence>)
    {
        val filePath = generateFilename(basePath, sample)

        val csvFormat = CSVFormat.Builder.create()
            .setDelimiter('\t').setRecordSeparator('\n')
            .setHeader(Column::class.java)
            .build()

        val sortedVdj = vdjSequences.sortedWith(
                Collections.reverseOrder(
                Comparator.comparingInt({ vdj: VDJSequence -> vdj.cdr3SupportMin })
                    .thenComparingInt({ vdj: VDJSequence -> vdj.numReads }) // handle the highest quality ones first
                    .thenComparingInt({ vdj: VDJSequence -> vdj.length })
                ))
        
        // in order to work out which ones are duplicate, we create a map of all VDJ sequences,
        // and the one with highest support count
        val cdr3SupportMap = HashMap<String, Int>()
        
        for (vdjSeq in sortedVdj)
        {
            val cdr3 = vdjSeq.cdr3Sequence
            cdr3SupportMap[cdr3] = Math.max(cdr3SupportMap.getOrDefault(cdr3, 0), vdjSeq.supportMin)
        }

        csvFormat.print(FileWriterUtils.createBufferedWriter(filePath)).use { printer: CSVPrinter ->
            for (vdj in sortedVdj)
            {
                val isDuplicate: Boolean = cdr3SupportMap.getOrDefault(vdj.cdr3Sequence, 0) > vdj.supportMin
                writeVDJSequence(printer, vdj, isDuplicate)
            }
        }
    }

    private fun writeVDJSequence(csvPrinter: CSVPrinter, vdj: VDJSequence, isDuplicate: Boolean)
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
                Column.vType -> csvPrinter.print(vdj.vAnchor.geneType)
                Column.vAnchorEnd -> csvPrinter.print(vdj.vAnchor.anchorBoundary)
                Column.vAnchorSeq -> csvPrinter.print(vdj.vAnchorSequence)
                Column.vAnchorTemplateSeq -> csvPrinter.print(vdj.vAnchor.templateAnchorSeq)
                Column.vAnchorAA -> csvPrinter.print(aminoAcidFromBases(vdj.vAnchorSequence))
                Column.vAnchorTemplateAA -> csvPrinter.print(aminoAcidFromBases(vdj.vAnchor.templateAnchorSeq))
                Column.vSimilarityScore -> csvPrinter.print(calcAnchorSimilarity(vdj, vdj.vAnchor))
                Column.jType -> csvPrinter.print(vdj.jAnchor.geneType)
                Column.jAnchorStart -> csvPrinter.print(vdj.jAnchor.anchorBoundary)
                Column.jAnchorSeq -> csvPrinter.print(vdj.jAnchorSequence)
                Column.jAnchorTemplateSeq -> csvPrinter.print(vdj.jAnchor.templateAnchorSeq)
                Column.jAnchorAA -> csvPrinter.print(aminoAcidFromBases(vdj.jAnchorSequence))
                Column.jAnchorTemplateAA -> csvPrinter.print(aminoAcidFromBases(vdj.jAnchor.templateAnchorSeq))
                Column.jSimilarityScore -> csvPrinter.print(calcAnchorSimilarity(vdj, vdj.jAnchor))
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

        if ((vdj.vAnchor is VJAnchorByBlosum) && (vdj.vAnchor.similarityScore < 0))
        {
            filters.add("NO_V_ANCHOR")
        }
        if ((vdj.jAnchor is VJAnchorByBlosum) && (vdj.jAnchor.similarityScore < 0))
        {
            filters.add("NO_J_ANCHOR")
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
}

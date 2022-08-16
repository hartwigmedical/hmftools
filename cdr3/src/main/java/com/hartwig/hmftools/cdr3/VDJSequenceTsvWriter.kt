package com.hartwig.hmftools.cdr3

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.FileWriterUtils
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import java.io.File
import java.util.*
import kotlin.collections.ArrayList

object VDJSequenceTsvWriter
{
    private enum class Column
    {
        id,
        filter,
        minHighQualBaseReads,
        assignedReads,
        inFrame,
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
        cdr3Seq,
        cdr3AA,
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

        csvFormat.print(FileWriterUtils.createBufferedWriter(filePath)).use { printer: CSVPrinter ->
            for (vdj in sortedVdj)
            {
                writeVDJSequence(printer, vdj)
            }
        }
    }

    private fun writeVDJSequence(csvPrinter: CSVPrinter, vdj: VDJSequence)
    {
        val vAnchorByBlosum: VJAnchorByBlosum? = vdj.vAnchor as? VJAnchorByBlosum
        val jAnchorByBlosum: VJAnchorByBlosum? = vdj.jAnchor as? VJAnchorByBlosum

        for (c in Column.values())
        {
            when (c)
            {
                Column.id -> csvPrinter.print(vdj.id)
                Column.filter -> csvPrinter.print(filterString(vdj).joinToString(separator = ","))
                Column.minHighQualBaseReads -> csvPrinter.print(vdj.cdr3SupportMin)
                Column.assignedReads -> csvPrinter.print(vdj.numReads)
                Column.inFrame -> csvPrinter.print((vdj.jAnchor.anchorBoundary % 3) == 0)
                Column.vType -> csvPrinter.print(vdj.vAnchor.geneType)
                Column.vAnchorEnd -> csvPrinter.print(vdj.vAnchor.anchorBoundary)
                Column.vAnchorSeq -> csvPrinter.print(vdj.vAnchorSequence)
                Column.vAnchorTemplateSeq -> csvPrinter.print(vAnchorByBlosum?.templateAnchorSeq ?: "null")
                Column.vAnchorAA -> csvPrinter.print(aminoAcidFromBases(vdj.vAnchorSequence))
                Column.vAnchorTemplateAA -> csvPrinter.print(if (vAnchorByBlosum != null)
                    aminoAcidFromBases(vAnchorByBlosum.templateAnchorSeq)
                    else "null")
                Column.vSimilarityScore -> csvPrinter.print(vAnchorByBlosum?.similarityScore ?: "null")
                Column.jType -> csvPrinter.print(vdj.jAnchor.geneType)
                Column.jAnchorStart -> csvPrinter.print(vdj.jAnchor.anchorBoundary)
                Column.jAnchorSeq -> csvPrinter.print(vdj.jAnchorSequence)
                Column.jAnchorTemplateSeq -> csvPrinter.print(jAnchorByBlosum?.templateAnchorSeq ?: "null")
                Column.jAnchorAA -> csvPrinter.print(aminoAcidFromBases(vdj.jAnchorSequence))
                Column.jAnchorTemplateAA -> csvPrinter.print(if (jAnchorByBlosum != null)
                    aminoAcidFromBases(jAnchorByBlosum.templateAnchorSeq)
                else "null")
                Column.jSimilarityScore -> csvPrinter.print(jAnchorByBlosum?.similarityScore ?: "null")
                Column.cdr3Seq -> csvPrinter.print(vdj.cdr3Sequence)
                Column.cdr3AA -> csvPrinter.print(aminoAcidFromBases(vdj.cdr3Sequence))
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

    fun filterString(vdj: VDJSequence): List<String>
    {
        val filters = ArrayList<String>()

        if ((vdj.jAnchor.anchorBoundary % 3) != 0)
        {
            filters.add("OUT_OF_FRAME")
        }
        if (vdj.aminoAcidSequence.contains(Codons.STOP_AMINO_ACID))
        {
            filters.add("CONTAINS_STOP")
        }
        if ((vdj.vAnchor is VJAnchorByBlosum) && (vdj.vAnchor.similarityScore < 0))
        {
            filters.add("NO_V_ANCHOR")
        }
        if ((vdj.jAnchor is VJAnchorByBlosum) && (vdj.jAnchor.similarityScore < 0))
        {
            filters.add("NO_J_ANCHOR")
        }

        if (filters.isEmpty())
            filters.add("PASS")

        return filters
    }
}

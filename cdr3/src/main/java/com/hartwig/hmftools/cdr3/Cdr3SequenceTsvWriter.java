package com.hartwig.hmftools.cdr3;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.cdr3.layout.VDJCandidate;
import com.hartwig.hmftools.common.codon.Codons;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.jetbrains.annotations.NotNull;

import kotlin.ranges.IntRange;

public class Cdr3SequenceTsvWriter
{
    private enum Column
    {
        FrameShifted, VGene, JGene, VAnchor, D, JAnchor, VAnchorAA, DAA, JAnchorAA, VLayoutSeq, JLayoutSeq, overlap
    }

    private static final String FILE_EXTENSION = ".cdr3.sequence.tsv";

    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(@NotNull final String filename,
            @NotNull final List<VDJCandidate> vdjCandidates) throws IOException
    {
        CSVFormat csvFormat = CSVFormat.Builder.create()
                .setDelimiter('\t').setRecordSeparator('\n')
                .setHeader(Column.class)
                .build();
        try (CSVPrinter csvPrinter = csvFormat.print(createBufferedWriter(filename)))
        {
            // for each key we need to find the V and the J and just put it down
            for (VDJCandidate vdj : vdjCandidates)
            {
                // we want to use the indices to work where things are
                String vdjSeq = vdj.getVdjSequence();
                List<VJGene> vGenes = VJReadLayoutAdaptor.getImgtGenes(vdj.getVLayout());
                IntRange vAnchorRange = VJReadLayoutAdaptor.getAnchorRange(VJGeneType.IGHV, vdj.getVLayout());

                List<VJGene> jGenes = VJReadLayoutAdaptor.getImgtGenes(vdj.getJLayout());
                // we need to shift it
                IntRange jAnchorRange = VJReadLayoutAdaptor.getAnchorRange(VJGeneType.IGHJ, vdj.getJLayout());
                int jSeqOffset = vdjSeq.length() - vdj.getJLayout().consensusSequence().length();
                jAnchorRange = new IntRange(jAnchorRange.getStart() + jSeqOffset, jAnchorRange.getEndInclusive() + jSeqOffset);

                String vAnchor = vdjSeq.substring(vAnchorRange.getStart(), vAnchorRange.getEndInclusive() + 1);
                String cdr = vdjSeq.substring(vAnchorRange.getEndInclusive() + 1, jAnchorRange.getStart());
                String jAnchor = vdjSeq.substring(jAnchorRange.getStart(), jAnchorRange.getEndInclusive() + 1);

                for (Column c : Column.values())
                {
                    switch (c)
                    {
                        case FrameShifted:
                            // NOTE: following are definitely a bit wrong, but we will sort it out later
                            csvPrinter.print((cdr.length() % 3) != 0);
                            break;
                        case VGene:
                            csvPrinter.print(vGenes.stream().map(VJGene::getName).distinct().collect(Collectors.toList()));
                            break;
                        case JGene:
                            csvPrinter.print(jGenes.stream().map(VJGene::getName).distinct().collect(Collectors.toList()));
                            break;
                        case VAnchor:
                            csvPrinter.print(vAnchor);
                            break;
                        case D:
                            csvPrinter.print(cdr);
                            break;
                        case JAnchor:
                            csvPrinter.print(jAnchor);
                            break;
                        case VAnchorAA:
                            csvPrinter.print(Codons.aminoAcidFromBases(vAnchor));
                            break;
                        case DAA:
                            csvPrinter.print(Codons.aminoAcidFromBases(cdr));
                            break;
                        case JAnchorAA:
                            csvPrinter.print(Codons.aminoAcidFromBases(jAnchor));
                            break;
                        case VLayoutSeq:
                            csvPrinter.print(vdj.getVLayout().consensusSequence());
                            break;
                        case JLayoutSeq:
                            csvPrinter.print(vdj.getJLayout().consensusSequence());
                            break;
                        case overlap:
                            csvPrinter.print(vdj.getOverlapSeq());
                            break;
                    }
                }
                csvPrinter.println();
            }
        }
    }
}

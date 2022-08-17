package com.hartwig.hmftools.cider;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;

import java.io.File;
import java.io.IOException;
import java.util.Collection;

import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.codon.Codons;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import org.jetbrains.annotations.NotNull;

public class Cdr3ReadTsvWriter
{
    private enum Column
    {
        readId, firstOfPair, chromosome, alignStart, alignEnd,
        vGene, jGene, vAnchorMatch, jAnchorMatch,
        vSoftClip, jSoftClip, dnaSeq, aaSeq
    }

    private static final String FILE_EXTENSION = ".cider.read.tsv";

    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    public static void write(@NotNull final String filename,
            @NotNull final Multimap<ReadKey, VJReadCandidate> cdr3Reads) throws IOException
    {
        CSVFormat csvFormat = CSVFormat.Builder.create()
                .setDelimiter('\t').setRecordSeparator('\n')
                .setHeader(Column.class)
                .build();
        try (CSVPrinter csvPrinter = csvFormat.print(createBufferedWriter(filename)))
        {
            // for each key we need to find the V and the J and just put it down
            for (ReadKey readKey : cdr3Reads.keySet())
            {
                Collection<VJReadCandidate> readCandidates = cdr3Reads.get(readKey);

                VJReadCandidate vMatch = readCandidates.stream()
                    .filter(o -> o.getVjGeneType() == VJGeneType.IGHV)
                    .findFirst()
                    .orElse(null);

                VJReadCandidate jMatch = readCandidates.stream()
                        .filter(o -> o.getVjGeneType() == VJGeneType.IGHJ)
                        .findFirst()
                        .orElse(null);

                if (vMatch == null && jMatch == null)
                    continue;

                Cdr3ReadVJMatch vjMatch = VJReadJoiner.combineVJMatch(readKey, vMatch, jMatch);

                // this is a bit wrong but we will live with it for now
                VJReadCandidate mainMatch = vjMatch.getMainMatch();

                csvPrinter.print(readKey.getReadName());
                csvPrinter.print(readKey.getFirstOfPair());
                csvPrinter.print(mainMatch.getRead().getReferenceName());
                csvPrinter.print(mainMatch.getRead().getAlignmentStart());
                csvPrinter.print(mainMatch.getRead().getAlignmentEnd());

                //csvPrinter.print(vMatch != null ? vMatch.getVjGenes().get(0).getName() : "none");
                //csvPrinter.print(jMatch != null ? jMatch.getVjGenes().get(0).getName() : "none");

                csvPrinter.print(vMatch != null ? vMatch.getAnchorMatchMethod() : "none");
                csvPrinter.print(jMatch != null ? jMatch.getAnchorMatchMethod() : "none");

                int vSoftClip = 0;
                int jSoftClip = 0;

                if (vMatch != null && vMatch.getAnchorMatchMethod() == VJReadCandidate.AnchorMatchMethod.ALIGN)
                {
                    // for v we want soft clip to be on the downstream side
                    vSoftClip = vMatch.getUseReverseComplement() ? vMatch.getLeftSoftClip() : vMatch.getRightSoftClip();
                }

                if (jMatch != null && jMatch.getAnchorMatchMethod() == VJReadCandidate.AnchorMatchMethod.ALIGN)
                {
                    // for v we want soft clip to be on the upstream side
                    jSoftClip = jMatch.getUseReverseComplement() ? jMatch.getRightSoftClip() : jMatch.getLeftSoftClip();
                }

                // NOTE: following are definitely a bit wrong, but we will sort it out later
                csvPrinter.print(vSoftClip);
                csvPrinter.print(jSoftClip);

                // now for the DNA sequence, we want to see if
                String cdr3Seq = vjMatch.getSequence();

                // ensure aligned to the codon
                int startIndex = mainMatch.getAnchorOffsetStart() % 3;
                int endIndex = cdr3Seq.length();

                if (vMatch != null)
                {
                    // we want the last 3 bases
                    startIndex = vjMatch.getVAnchorOffsetEnd() - 3;
                }

                if (jMatch != null)
                {
                    endIndex = vjMatch.getJAnchorOffsetStart() + 3;
                }

                cdr3Seq = cdr3Seq.substring(startIndex, endIndex);

                csvPrinter.print(cdr3Seq);
                csvPrinter.print(Codons.aminoAcidFromBases(cdr3Seq));
                csvPrinter.println();
            }
        }
    }
}

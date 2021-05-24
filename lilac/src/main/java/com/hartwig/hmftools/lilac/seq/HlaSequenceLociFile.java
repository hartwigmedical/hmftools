package com.hartwig.hmftools.lilac.seq;

import static java.lang.Math.max;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;

public class HlaSequenceLociFile
{
    // writes exon boundaries, and then the allele sequences
    public static void write(
            final String fileName, final List<Integer> aBoundaries, final List<Integer> bBoundaries,
            final List<Integer> cBoundaries, final List<HlaSequenceLoci> sequences)
    {
        /*
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.newLine();



            writer.close();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", fileName, e.toString());
            return;
        }

        val outputFile = File(file)
        outputFile.writeText("")

        val maxLengths = sequences.map { it.lengths() }.reduce { left, right -> maxLengths(left, right) }
        val templateSequence = sequences[0].padInserts(maxLengths)
        outputFile.appendText("HLA-A Boundaries".padEnd(20, ' ') + "\t" + boundaryString(aBoundaries, maxLengths) + "\n")
        outputFile.appendText("HLA-B Boundaries".padEnd(20, ' ') + "\t" + boundaryString(bBoundaries, maxLengths) + "\n")
        outputFile.appendText("HLA-C Boundaries".padEnd(20, ' ') + "\t" + boundaryString(cBoundaries, maxLengths) + "\n")

        outputFile.appendText("${sequences[0].allele}".padEnd(20, ' ') + "\t" + templateSequence + "\n")
        for (i in 1 until sequences.size) {
            val victimSequence = diff(sequences[i].padInserts(maxLengths), templateSequence)
            outputFile.appendText("${sequences[i].allele}".padEnd(20, ' ') + "\t" + victimSequence + "\n")
        }
         */
    }

    private final String boundaryString(Collection<Integer> boundaries, List<Integer> lengths)
    {
        /*
        val range = (0..boundaries.max()!!)
        val unpadded = range.mapIndexed { index, _ -> if (index in boundaries) "|" else " " }
        return unpadded.mapIndexed { index, value -> value.padEnd(lengths[index], ' ') }.joinToString("")

         */
        return "";
    }

    private final String diff(String victim, String reference)
    {
        /*
        return victim
        .mapIndexed { index, _ -> diff(victim[index], if (index < reference.length) reference[index] else '!') }
        .joinToString("")

         */
        return "";
    }

    private final char diff(char c, char reference)
    {
        return (c == '|' ? 124 : (c == '*' ? 42 : (c == '.' ? 46 : (c == reference ? 45 : c))));
    }

    private String padInserts(final HlaSequenceLoci seqLoci, final List<Integer> lengths)
    {
        for(int i = 0; i < seqLoci.getSequences().size(); ++i)
        {
            int reqLength = lengths.get(i);
            String insert = seqLoci.getSequences().get(i);

            for(int j = insert.length(); j < lengths.get(i); ++j)
            {
                insert += ".";
            }
        }

        /*
        return this.sequences
                .mapIndexed { index, value -> value.padEnd(lengths[index], '.') }
                .joinToString("")
            .trimEnd { x -> x == '.' }
         */

        return "";
    }

    private List<Integer> lengths(final HlaSequenceLoci seqLoci)
    {
        return seqLoci.getSequences().stream().map(x -> new Integer(x.length())).collect(Collectors.toList());
    }

    private List<Integer> maxLengths(final List<Integer> left, final List<Integer> right)
    {
        List<Integer> results = Lists.newArrayList();

        for (int i = 0; i <  max(left.size(), right.size()); ++i)
        {
            int leftValue = i < left.size() ? left.get(i) : 0;
            int rightValue = i < right.size() ? right.get(i) : 0;
            results.add(max(leftValue, rightValue));
        }

        return results;
    }

}

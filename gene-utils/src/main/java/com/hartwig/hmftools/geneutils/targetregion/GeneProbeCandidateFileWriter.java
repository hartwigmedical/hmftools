package com.hartwig.hmftools.geneutils.targetregion;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

import org.apache.commons.lang3.tuple.Pair;

public class GeneProbeCandidateFileWriter
{
    enum Column
    {
        GeneName,
        RegionType,
        Chromosome,
        RegionStart,
        RegionEnd,
        ProbeStart,
        ProbeEnd,
        ProbeGcContent,
        ProbeSumBitScore,
        Selected,
        ProbeSequence
    }

    private static final String FILE_EXTENSION = ".probe_candidate.tsv";

    public static String generateFilename(String basePath, String prefix)
    {
        return basePath + File.separator + prefix + FILE_EXTENSION;
    }

    public static void write(final String basePath, final String prefix, final List<TargetedGeneRegion> targetedGeneRegions)
    {
        final String fileName = generateFilename(basePath, prefix);

        List<Pair<TargetedGeneRegion, ProbeCandidate>> probeList = targetedGeneRegions.stream()
                .flatMap(o -> o.getProbeCandidates().stream().map(c -> Pair.of(o, c)))
                .collect(Collectors.toList());

        try(BufferedWriter writer = createBufferedWriter(fileName))
        {
            DelimFileWriter.write(writer, Column.values(), probeList,
                    (r, row) -> {
                        row.set(Column.GeneName, r.getLeft().getGene().getGeneData().GeneName);
                        row.set(Column.RegionType, r.getLeft().getType().name());
                        row.set(Column.Chromosome, r.getLeft().getChromosome());
                        row.set(Column.RegionStart, r.getLeft().getStart());
                        row.set(Column.RegionEnd, r.getLeft().getEnd());
                        row.set(Column.ProbeStart, r.getRight().getStart());
                        row.set(Column.ProbeEnd, r.getRight().getEnd());
                        row.set(Column.ProbeGcContent, r.getRight().getGcContent());
                        row.set(Column.ProbeSumBitScore, r.getRight().getSumBlastnBitScore());
                        row.set(Column.Selected, Boolean.toString(r.getLeft().getSelectedProbe() == r.getRight()));
                        row.set(Column.ProbeSequence, r.getRight().getSequence());
                    });
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }
}

package com.hartwig.hmftools.geneutils.targetregion;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.List;

import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

public class GeneProbeRegionFileWriter
{
    enum Column
    {
        GeneName,
        RegionType,
        Chromosome,
        RegionStart,
        RegionEnd,
        UseWholeRegion,
        ProbeStart,
        ProbeEnd,
        ProbeGcContent,
        ProbeSumBitScore,
        ProbeSequence
    }

    private static final String FILE_EXTENSION = ".gene_region.tsv";

    public static String generateFilename(String basePath, String prefix)
    {
        return basePath + File.separator + prefix + FILE_EXTENSION;
    }

    public static void write(final String basePath, final String prefix, final List<TargetedGeneRegion> targetedGeneRegions)
    {
        final String fileName = generateFilename(basePath, prefix);

        try(BufferedWriter writer = createBufferedWriter(fileName))
        {
            DelimFileWriter.write(writer, Column.values(), targetedGeneRegions,
                (r, row) -> {
                    row.set(Column.GeneName, r.getGene().getGeneData().GeneName);
                    row.set(Column.RegionType, r.getType().name());
                    row.set(Column.Chromosome, r.getChromosome());
                    row.set(Column.RegionStart, r.getStart());
                    row.set(Column.RegionEnd, r.getEnd());
                    row.set(Column.UseWholeRegion, r.useWholeRegion());

                    // some we use whole region
                    if(!r.useWholeRegion())
                    {
                        ProbeCandidate selectedProbe = r.getSelectedProbe();
                        if(selectedProbe != null)
                        {
                            row.set(Column.ProbeStart, selectedProbe.getStart());
                            row.set(Column.ProbeEnd, selectedProbe.getEnd());
                            row.set(Column.ProbeGcContent, selectedProbe.getGcContent());
                            row.set(Column.ProbeSumBitScore, selectedProbe.getSumBlastnBitScore());
                        }
                    }
                });
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }
}

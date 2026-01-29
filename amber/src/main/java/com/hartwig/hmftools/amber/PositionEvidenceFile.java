package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createGzipBufferedWriter;

import java.io.IOException;
import java.io.Writer;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.ListMultimap;
import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

public final class PositionEvidenceFile
{

    private static final String EXTENSION = ".amber.raw.tsv.gz";

    public static String generateAmberFilenameForWriting(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    private static final String CHROMOSOME = "chromosome";
    private static final String POSITION = "position";
    private static final String REF = "ref";
    private static final String ALT = "alt";
    private static final String READ_DEPTH = "readDepth";
    private static final String INDEL_COUNT = "indelCount";
    private static final String REF_COUNT = "refCount";
    private static final String ALT_COUNT = "altCount";

    public static void write(String filename, ListMultimap<Chromosome, PositionEvidence> data) throws IOException
    {
        List<PositionEvidence> bafs = new ArrayList<>();
        for(Chromosome chromosome : HumanChromosome.values())
        {
            if(chromosome != HumanChromosome._11)
            {
                continue;
            }

            if(data.containsKey(chromosome))
            {
                bafs.addAll(data.get(chromosome).stream().filter(p -> p.position() < 10_000_001).sorted().toList());
            }
        }
        write(filename, bafs);
    }

    public static void write(final String filename, final List<PositionEvidence> bafs) throws IOException
    {
        try(Writer writer = createGzipBufferedWriter(filename))
        {
            for(String line : toLines(bafs))
            {
                writer.write(line + '\n');
            }
        }
    }

    private static List<String> toLines(final List<PositionEvidence> bafs)
    {
        final List<String> lines = new ArrayList<>();
        lines.add(header());
        bafs.stream().map(PositionEvidenceFile::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM, "", "")
                .add(CHROMOSOME)
                .add(POSITION)
                .add(REF)
                .add(ALT)
                .add(READ_DEPTH)
                .add(INDEL_COUNT)
                .add(REF_COUNT)
                .add(ALT_COUNT)
                .toString();
    }

    private static String toString(final PositionEvidence baf)
    {
        return new StringJoiner(TSV_DELIM)
                .add(baf.Chromosome)
                .add(String.valueOf(baf.Position))
                .add(baf.ref())
                .add(baf.alt())
                .add(String.valueOf(baf.ReadDepth))
                .add(String.valueOf(baf.IndelCount))
                .add(String.valueOf(baf.RefSupport))
                .add(String.valueOf(baf.AltSupport))
                .toString();
    }
}

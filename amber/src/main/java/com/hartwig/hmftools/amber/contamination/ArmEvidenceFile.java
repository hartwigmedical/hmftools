package com.hartwig.hmftools.amber.contamination;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.segmentation.ChrArm;

public final class ArmEvidenceFile
{
    private static final String EXTENSION = ".amber.contamination.arm.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + EXTENSION;
    }

    private static final String CHROMOSOME = "chromosome";
    private static final String ARM = "arm";
    private static final String TOTAL_POINTS = "totalPoints";
    private static final String EVIDENCE_POINTS = "evidencePoints";

    public static void write(final String filename, final Collection<CategoryEvidence<ChrArm>> categoryEvidences) throws IOException
    {
        try(Writer writer = createBufferedWriter(filename))
        {
            for(String line : toLines(categoryEvidences))
            {
                writer.write(line + '\n');
            }
        }
    }

    private static List<String> toLines(final Collection<CategoryEvidence<ChrArm>> categoryEvidences)
    {
        final List<String> lines = new ArrayList<>();
        lines.add(header());
        categoryEvidences.stream().sorted().map(ArmEvidenceFile::toString).forEach(lines::add);
        return lines;
    }

    private static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add(CHROMOSOME)
                .add(ARM)
                .add(EVIDENCE_POINTS)
                .add(TOTAL_POINTS)
                .toString();
    }

    private static String toString(final CategoryEvidence<ChrArm> categoryEvidence)
    {
        return new StringJoiner(TSV_DELIM)
                .add(categoryEvidence.category().chromosome().toString())
                .add(categoryEvidence.category().arm().toString())
                .add(String.valueOf(categoryEvidence.evidencePoints()))
                .add(String.valueOf(categoryEvidence.totalPoints()))
                .toString();
    }
}

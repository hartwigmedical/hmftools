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

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.segmentation.Arm;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.common.utils.file.DelimFileWriter;

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

    private static List<String> columns()
    {
        return List.of(CHROMOSOME, ARM, TOTAL_POINTS, EVIDENCE_POINTS);
    }

    public static void write(final String filename, final Collection<CategoryEvidence<ChrArm>> categoryEvidences) throws IOException
    {
        List<CategoryEvidence<ChrArm>> sorted = categoryEvidences.stream().sorted().toList();
        DelimFileWriter.write(filename, columns(), sorted, (baf, row) ->
        {
            row.set(CHROMOSOME, baf.category().chromosome().shortName());
            row.set(ARM, baf.category().arm().name());
            row.set(TOTAL_POINTS, baf.totalPoints());
            row.set(EVIDENCE_POINTS, baf.evidencePoints());
        });
    }

    public static List<CategoryEvidence<ChrArm>> read(String filename)
    {
        List<CategoryEvidence<ChrArm>> result = new ArrayList<>();
        try(DelimFileReader reader = new DelimFileReader(filename))
        {
            for(DelimFileReader.Row row : reader)
            {
                HumanChromosome chromosome = HumanChromosome.fromString(row.get(CHROMOSOME));
                Arm arm = Arm.fromString(row.get(ARM));
                int totalPoints = Integer.parseInt(row.get(TOTAL_POINTS));
                int evidencePoints = Integer.parseInt(row.get(EVIDENCE_POINTS));
                CategoryEvidence<ChrArm> evidence = new CategoryEvidence<ChrArm>(new ChrArm(chromosome, arm));
                evidence.set(totalPoints, evidencePoints);
                result.add(evidence);
            }
        }
        return result;
    }
}

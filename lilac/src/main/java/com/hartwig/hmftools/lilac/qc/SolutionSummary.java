package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.DELIM;

import com.hartwig.hmftools.lilac.coverage.AlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverage;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

public class SolutionSummary
{
    public final ComplexCoverage ReferenceCoverage;
    public final ComplexCoverage TumorCoverage;
    public final List<Double> TumorCopyNumber;
    public final List<SomaticCodingCount> SomaticCodingCount;
    public final ComplexCoverage RnaCoverage;

    public SolutionSummary(
            final ComplexCoverage referenceCoverage, final ComplexCoverage tumorCoverage,
            final List<Double> tumorCopyNumber, final List<SomaticCodingCount> somaticCodingCount, final ComplexCoverage rnaCoverage)
    {
        ReferenceCoverage = referenceCoverage;
        TumorCoverage = tumorCoverage;
        TumorCopyNumber = tumorCopyNumber;
        SomaticCodingCount = somaticCodingCount;
        RnaCoverage = rnaCoverage;
    }

    public final void write(final String fileName)
    {
        try
        {
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write(generateAlleleHeader());
            writer.newLine();

            for(int i = 0; i < 6; ++i)
            {
                writer.write(generateAlleleBody(i));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", fileName, e.toString());
            return;
        }
    }

    private final String generateAlleleHeader()
    {
        StringJoiner header = new StringJoiner(DELIM)
                .add("Allele")
                .add("RefTotal")
                .add("RefUnique")
                .add("RefShared")
                .add("RefWild")
                .add("TumorTotal")
                .add("TumorUnique")
                .add("TumorShared")
                .add("TumorWild")
                .add("TumorCopyNumber")
                .add("RnaTotal")
                .add("RnaUnique")
                .add("RnaShared")
                .add("RnaWild")
                .add("SomaticMissense")
                .add("SomaticNonsenseOrFrameshift")
                .add("SomaticSplice")
                .add("SomaticSynonymous")
                .add("SomaticInframeIndel");
        return header.toString();
    }

    private final String generateAlleleBody(int index)
    {
        AlleleCoverage ref = ReferenceCoverage.getAlleleCoverage().get(index);

        AlleleCoverage tumor = !TumorCoverage.getAlleleCoverage().isEmpty() ?
                TumorCoverage.getAlleleCoverage().get(index) : new AlleleCoverage(ref.Allele, 0, 0, 0);

        AlleleCoverage rna = !RnaCoverage.getAlleleCoverage().isEmpty() ?
                RnaCoverage.getAlleleCoverage().get(index) : new AlleleCoverage(ref.Allele, 0, 0, 0);

        double copyNumber = TumorCopyNumber.get(index);
        SomaticCodingCount codingCount = SomaticCodingCount.get(index);

        StringJoiner header = new StringJoiner(DELIM)
                .add(ref.Allele.toString())
                .add(String.valueOf(round(ref.TotalCoverage)))
                .add(String.valueOf(ref.UniqueCoverage))
                .add(String.valueOf(round(ref.SharedCoverage)))
                .add(String.valueOf(round(ref.WildCoverage)))
                .add(String.valueOf(round(tumor.TotalCoverage)))
                .add(String.valueOf(tumor.UniqueCoverage))
                .add(String.valueOf(round(tumor.SharedCoverage)))
                .add(String.valueOf(round(tumor.WildCoverage)))
                .add(String.format("%.2f", copyNumber))
                .add(String.valueOf(round(rna.TotalCoverage)))
                .add(String.valueOf(rna.UniqueCoverage))
                .add(String.valueOf(round(rna.SharedCoverage)))
                .add(String.valueOf(round(rna.WildCoverage)))
                .add(String.format("%.2f", codingCount.missense()))
                .add(String.format("%.2f", codingCount.nonsense()))
                .add(String.format("%.2f", codingCount.splice()))
                .add(String.format("%.2f", codingCount.synonymous()))
                .add(String.format("%.2f", codingCount.inframeIndel()));

        return header.toString();
    }

    public static SolutionSummary create(
            final ComplexCoverage referenceCoverage, final ComplexCoverage tumorCoverage,
            final List<Double> tumorCopyNumber, final List<SomaticCodingCount> somaticCodingCount, final ComplexCoverage rnaCoverage)
    {
        List<SomaticCodingCount> sortedCodingCount = somaticCodingCount.stream().collect(Collectors.toList());
        Collections.sort(sortedCodingCount, new SomaticCodingCountSorter());

        return new SolutionSummary(referenceCoverage, tumorCoverage, tumorCopyNumber, sortedCodingCount, rnaCoverage);
    }

    private static class SomaticCodingCountSorter implements Comparator<SomaticCodingCount>
    {
        public int compare(final SomaticCodingCount first, final SomaticCodingCount second)
        {
            return first.Allele.compareTo(second.Allele);
        }
    }

}

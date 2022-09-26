package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.round;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;

import com.hartwig.hmftools.common.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.lilac.coverage.AlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverage;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.compress.utils.Lists;

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
            List<LilacAllele> alleles = Lists.newArrayList();

            if(ReferenceCoverage != null)
            {
                for(int i = 0; i < 6; ++i)
                {
                    alleles.add(buildAlleleData(i));
                }
            }

            LilacAllele.write(fileName, alleles);
        }
        catch(IOException e)
        {
            LL_LOGGER.error("failed to write {}: {}", fileName, e.toString());
            return;
        }
    }

    private LilacAllele buildAlleleData(int index)
    {
        AlleleCoverage ref = ReferenceCoverage.getAlleleCoverage().get(index);

        AlleleCoverage tumor = !TumorCoverage.getAlleleCoverage().isEmpty() ?
                TumorCoverage.getAlleleCoverage().get(index) : new AlleleCoverage(ref.Allele, 0, 0, 0);

        AlleleCoverage rna = !RnaCoverage.getAlleleCoverage().isEmpty() ?
                RnaCoverage.getAlleleCoverage().get(index) : new AlleleCoverage(ref.Allele, 0, 0, 0);

        double copyNumber = TumorCopyNumber.get(index);
        SomaticCodingCount codingCount = SomaticCodingCount.get(index);

        return ImmutableLilacAllele.builder()
                .allele(ref.Allele.toString())
                .refFragments((int)round(ref.TotalCoverage))
                .refUnique(ref.UniqueCoverage)
                .refShared((int)round(ref.SharedCoverage))
                .refWild((int)round(ref.WildCoverage))
                .tumorFragments((int)round(tumor.TotalCoverage))
                .tumorUnique(tumor.UniqueCoverage)
                .tumorShared((int)round(tumor.SharedCoverage))
                .tumorWild((int)round(tumor.WildCoverage))
                .tumorCopyNumber(copyNumber)
                .rnaFragments((int)round(rna.TotalCoverage))
                .rnaUnique(rna.UniqueCoverage)
                .rnaShared((int)round(rna.SharedCoverage))
                .rnaWild((int)round(rna.WildCoverage))
                .somaticMissense(codingCount.missense())
                .somaticNonsenseOrFrameshift(codingCount.nonsense())
                .somaticSplice(codingCount.splice())
                .somaticSynonymous(codingCount.synonymous())
                .somaticInframeIndel(codingCount.inframeIndel())
                .build();
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

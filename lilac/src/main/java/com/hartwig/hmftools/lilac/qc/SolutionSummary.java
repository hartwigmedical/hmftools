package com.hartwig.hmftools.lilac.qc;

import static java.lang.Math.round;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.CURRENT_GENES;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.utils.file.FileLock;
import com.hartwig.hmftools.lilac.GeneSelector;
import com.hartwig.hmftools.lilac.coverage.AlleleCoverage;
import com.hartwig.hmftools.lilac.coverage.ComplexCoverage;
import com.hartwig.hmftools.lilac.hla.HlaAllele;
import com.hartwig.hmftools.lilac.hla.HlaGene;
import com.hartwig.hmftools.lilac.variant.SomaticCodingCount;

public class SolutionSummary
{
    private final GeneSelector mGenes;

    public final ComplexCoverage ReferenceCoverage;
    public final ComplexCoverage TumorCoverage;
    public final List<Double> TumorCopyNumber;
    public final List<SomaticCodingCount> SomaticCodingCount;
    public final ComplexCoverage RnaCoverage;

    public SolutionSummary(
            final GeneSelector genes, final ComplexCoverage referenceCoverage, final ComplexCoverage tumorCoverage,
            final List<Double> tumorCopyNumber, final List<SomaticCodingCount> somaticCodingCount, final ComplexCoverage rnaCoverage)
    {
        mGenes = genes;
        ReferenceCoverage = referenceCoverage;
        TumorCoverage = tumorCoverage;
        TumorCopyNumber = tumorCopyNumber;
        SomaticCodingCount = somaticCodingCount;
        RnaCoverage = rnaCoverage;
    }

    private LilacAllele buildAlleleData(int index)
    {
        // ref will be empty in tumor-only mode
        HlaAllele refAllele = !ReferenceCoverage.getAlleleCoverage().isEmpty() ?
                ReferenceCoverage.getAlleleCoverage().get(index).Allele : TumorCoverage.getAlleleCoverage().get(index).Allele;

        AlleleCoverage noCoverage = new AlleleCoverage(refAllele, 0, 0, 0);

        AlleleCoverage ref = !ReferenceCoverage.getAlleleCoverage().isEmpty() ? ReferenceCoverage.getAlleleCoverage().get(index) : noCoverage;

        AlleleCoverage tumor = !TumorCoverage.getAlleleCoverage().isEmpty() ? TumorCoverage.getAlleleCoverage().get(index) : noCoverage;

        AlleleCoverage rna = !RnaCoverage.getAlleleCoverage().isEmpty() ? RnaCoverage.getAlleleCoverage().get(index) : noCoverage;

        double copyNumber = TumorCopyNumber.get(index);
        SomaticCodingCount codingCount = SomaticCodingCount.get(index);

        return ImmutableLilacAllele.builder()
                .allele(refAllele.toString())
                .refFragments((int) round(ref.TotalCoverage))
                .refUnique(ref.UniqueCoverage)
                .refShared((int) round(ref.SharedCoverage))
                .refWild((int) round(ref.WildCoverage))
                .tumorFragments((int) round(tumor.TotalCoverage))
                .tumorUnique(tumor.UniqueCoverage)
                .tumorShared((int) round(tumor.SharedCoverage))
                .tumorWild((int) round(tumor.WildCoverage))
                .tumorCopyNumber(copyNumber)
                .rnaFragments((int) round(rna.TotalCoverage))
                .rnaUnique(rna.UniqueCoverage)
                .rnaShared((int) round(rna.SharedCoverage))
                .rnaWild((int) round(rna.WildCoverage))
                .somaticMissense(codingCount.missense())
                .somaticNonsenseOrFrameshift(codingCount.nonsense())
                .somaticSplice(codingCount.splice())
                .somaticSynonymous(codingCount.synonymous())
                .somaticInframeIndel(codingCount.inframeIndel())
                .build();
    }

    public static SolutionSummary create(
            final ComplexCoverage referenceCoverage, final ComplexCoverage tumorCoverage,
            final List<Double> tumorCopyNumber, final Iterable<SomaticCodingCount> somaticCodingCount, final ComplexCoverage rnaCoverage)
    {
        List<SomaticCodingCount> sortedCodingCount = Lists.newArrayList(somaticCodingCount);
        Collections.sort(sortedCodingCount, new SomaticCodingCountSorter());

        return new SolutionSummary(CURRENT_GENES, referenceCoverage, tumorCoverage, tumorCopyNumber, sortedCodingCount, rnaCoverage);
    }

    public void write(final String fileName)
    {
        List<LilacAllele> alleles = Lists.newArrayList();
        Set<HlaGene> genes = Sets.newHashSet();
        if(ReferenceCoverage != null)
        {
            for(int i = 0; i < ReferenceCoverage.getAlleles().size(); ++i)
            {
                genes.add(ReferenceCoverage.getAlleles().get(i).Gene);
                alleles.add(buildAlleleData(i));
            }
        }

        File file = new File(fileName);
        List<String> existingLines = Lists.newArrayList();
        try(FileLock fileLock = FileLock.create(file))
        {
            BufferedReader reader = fileLock.getBufferedReader();
            reader.readLine();
            String line = reader.readLine();
            while(line != null)
            {
                String geneStr = line.split("\\*")[0];
                HlaGene gene = HlaGene.fromString(geneStr);
                if(gene != null && !mGenes.contains(gene))
                    existingLines.add(line);

                line = reader.readLine();
            }

            fileLock.clear();
            BufferedWriter writer = fileLock.getBufferedWriter();
            LilacAllele.write(writer, alleles, existingLines);
            writer.flush();
        }
        catch(Exception e)
        {
            LL_LOGGER.error("failed to update {}: {}", fileName, e.toString());
        }
    }

    private static class SomaticCodingCountSorter implements Comparator<SomaticCodingCount>
    {
        public int compare(final SomaticCodingCount first, final SomaticCodingCount second)
        {
            return first.Allele.compareTo(second.Allele);
        }
    }
}
